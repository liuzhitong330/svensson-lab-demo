"""Exploratory precursor-RNA prioritization; no inference of mature peptide abundance.

Run fetch_inputs.py first. BioMart mapping downloads are specified in the README.
Dependencies: numpy, pandas, openpyxl, regex. All upstream code stays in private cache.
"""
from pathlib import Path
from types import SimpleNamespace
from collections import defaultdict, Counter
import argparse
import ast
import concurrent.futures
import hashlib
import json
import math
import re
import sys
import time
import urllib.request
import numpy as np
import pandas as pd

GROUPS = {"LFD": ["B1", "B2", "B4", "B5"], "HFD": ["D1", "D2", "D3", "D5"],
          "SEMA": ["F2", "F4", "F5"]}
SAMPLES = sum(GROUPS.values(), [])
MAIN = GROUPS["HFD"] + GROUPS["SEMA"]
EXPECTED = {"FastaParser.py": "1ad36e296cf66e7d488fed5cdd765cb1c8ddb07bb27d4c8abff66527e4cf5436",
            "secreted.fasta": "35cad8ba94f0a201f706826a789adabfc744023132f265d368b27988875be6d1"}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ratio_factors(frame):
    """DESeq-style median ratios on genes strictly positive in every included sample."""
    positive = frame.loc[(frame > 0).all(axis=1)].astype(float)
    if positive.empty:
        raise ValueError("No all-positive genes for median-ratio normalization")
    geometric = np.exp(np.log(positive).mean(axis=1))
    factors = positive.div(geometric, axis=0).median(axis=0)
    assert np.isfinite(factors).all() and (factors > 0).all()
    return factors, len(positive)


def effect(row, groups=GROUPS):
    return float(np.log2((row[groups["SEMA"]].mean()+1)/(row[groups["HFD"]].mean()+1)))


def parse_fasta(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            accession = line.split("|")[1]
            result[accession] = {"sequence": "", "header": line[1:],
                                 "gene": re.search(r"GN=(\S+)", line).group(1) if "GN=" in line else None}
        else:
            result[accession]["sequence"] += line.strip()
    return result


def predictor(cache):
    """Run only the hash-verified upstream expression + pure cleavage finder."""
    for name, expected in EXPECTED.items():
        if digest(cache/name) != expected:
            raise ValueError(f"Unreviewed upstream version: {name}")
    try:
        import regex
    except ImportError:
        sys.path.insert(0, str(cache.parent/"python-deps"))
        import regex
    tree = ast.parse((cache/"FastaParser.py").read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef))
    init = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "__init__")
    assignment = next(n for n in init.body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Attribute) and t.attr == "prohormone_regex" for t in n.targets))
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "get_cleavage_sites")
    namespace = {"regex": regex}
    exec(compile(ast.Module(body=[method], type_ignores=[]), "reviewed-upstream-method", "exec"), namespace)
    instance = SimpleNamespace(cleavage_spacing=5)
    instance.prohormone_regex = eval(compile(ast.Expression(assignment.value), "reviewed-upstream-expression", "eval"),
                                    {"constants": SimpleNamespace(SIGNAL_PEPTIDE_LENGTH=20), "self": instance})
    def cuts(sequence):
        matches = namespace["get_cleavage_sites"](instance, SimpleNamespace(sequence_as_string=lambda: sequence))
        return {(m.start(), m.end()) for m in matches}
    return cuts


def fetch_homology(cache, gene):
    path = cache/"homology"/f"{gene}.json"
    if path.exists():
        return json.loads(path.read_text())
    url = f"https://rest.ensembl.org/homology/id/human/{gene}?target_species=mus_musculus;type=orthologues;content-type=application/json"
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "ResearchReanalysis/1.0"})
            with urllib.request.urlopen(request, timeout=45) as response:
                data = json.load(response)
            path.write_text(json.dumps(data))
            return data
        except Exception as exc:
            if attempt == 3:
                return {"fetchError": str(exc), "url": url}
            time.sleep(1.5*(attempt+1))


def peptide_projection(peptide, human, source_alignment, target_alignment, cuts):
    """Project exact, uniquely placed human peptide through actual Ensembl alignment."""
    sequence = peptide["sequence"]
    starts = [m.start() for m in re.finditer(f"(?={re.escape(sequence)})", human)]
    base = {"humanSequence": sequence, "sourceRow": peptide["sourceRow"],
            "reportedLength": peptide["reportedLength"], "actualLength": len(sequence),
            "sequenceConserved": False, "cleavageConserved": False}
    if len(starts) != 1:
        return {**base, "status": "ambiguous_or_missing_human_location"}
    start, end = starts[0], starts[0]+len(sequence)
    if start < 2 or end+2 > len(human):
        return {**base, "status": "terminal_fragment_without_two_flanks"}
    source_columns = [i for i, aa in enumerate(source_alignment) if aa != "-"]
    col_start, col_end = source_columns[start], source_columns[end-1]+1
    aligned_h = source_alignment[col_start:col_end]
    aligned_m = target_alignment[col_start:col_end]
    mouse_sequence = target_alignment.replace("-", "")
    mouse_start = len(target_alignment[:col_start].replace("-", ""))
    mouse_end = len(target_alignment[:col_end].replace("-", ""))
    mouse_peptide = target_alignment[col_start:col_end].replace("-", "")
    human_sites, mouse_sites = cuts(human), cuts(mouse_sequence)
    human_boundaries = ((start-2, start) in human_sites and (end, end+2) in human_sites
                        and not any(start <= a and b <= end for a, b in human_sites))
    mouse_boundaries = ((mouse_start-2, mouse_start) in mouse_sites and (mouse_end, mouse_end+2) in mouse_sites
                        and not any(mouse_start <= a and b <= mouse_end for a, b in mouse_sites))
    return {
        **base, "status": "projected", "humanStart": start+1, "humanEnd": end,
        "mouseStart": mouse_start+1, "mouseEnd": mouse_end,
        "mouseSequence": mouse_peptide, "alignedHuman": aligned_h, "alignedMouse": aligned_m,
        "alignmentHasIndel": "-" in aligned_h or "-" in aligned_m,
        "leftFlankHuman": human[start-2:start], "rightFlankHuman": human[end:end+2],
        "leftFlankMouse": mouse_sequence[max(0, mouse_start-2):mouse_start],
        "rightFlankMouse": mouse_sequence[mouse_end:mouse_end+2],
        "sequenceConserved": sequence == mouse_peptide and "-" not in aligned_h and "-" not in aligned_m,
        "cleavageConserved": human_boundaries and mouse_boundaries,
        "humanBoundariesPredicted": human_boundaries, "mouseBoundariesPredicted": mouse_boundaries,
        "identicalResidues": sum(a == b and a != "-" for a, b in zip(aligned_h, aligned_m)),
        "alignmentLength": len(aligned_h),
        "mismatchPositions": [i+1 for i, (a, b) in enumerate(zip(aligned_h, aligned_m)) if a != b],
    }


def main(cache, output):
    cache.mkdir(parents=True, exist_ok=True)
    (cache/"homology").mkdir(exist_ok=True)
    counts = pd.read_csv(cache/"GSE266899_raw_counts.txt.gz", sep="\t", index_col=0)
    counts.index = counts.index.str.replace(r"\.\d+$", "", regex=True)
    assert counts.index.is_unique, "Version stripping must not silently collapse rows"
    assert (counts >= 0).all().all() and np.equal(counts, np.floor(counts)).all().all()
    t2 = counts[SAMPLES]
    sf, normalization_genes = ratio_factors(t2)
    normalized = t2.div(sf)
    cpm = t2.div(t2.sum())*1e6
    effects = np.log2((normalized[GROUPS["SEMA"]].mean(axis=1)+1)/(normalized[GROUPS["HFD"]].mean(axis=1)+1))
    cpm_effects = np.log2((cpm[GROUPS["SEMA"]].mean(axis=1)+1)/(cpm[GROUPS["HFD"]].mean(axis=1)+1))
    loo = {}
    for omit in MAIN:
        reduced = t2.drop(columns=omit)
        reduced_sf, _ = ratio_factors(reduced)
        reduced_norm = reduced.div(reduced_sf)
        g = {name: [s for s in samples if s != omit] for name, samples in GROUPS.items()}
        loo[omit] = np.log2((reduced_norm[g["SEMA"]].mean(axis=1)+1)/(reduced_norm[g["HFD"]].mean(axis=1)+1))
    print("Normalization complete", normalization_genes, sf.to_dict(), flush=True)
    fasta = parse_fasta(cache/"secreted.fasta")
    published = json.loads((cache/"published-peptides.json").read_text())
    short = [p for p in published if 5 <= p["actualLength"] <= 25]
    by_protein = defaultdict(list)
    for p in short:
        by_protein[p["uniprot"]].append(p)
    ortholog = pd.read_csv(cache/"human-mouse-biomart.tsv", sep="\t", dtype=str).fillna("")
    mapping_path = cache/"human-identifiers-biomart.tsv"
    if not mapping_path.exists():
        alternatives = [p for p in cache.glob("*biomart*.tsv") if "hgnc_symbol" in p.read_text().splitlines()[0]]
        if not alternatives:
            raise RuntimeError("Awaiting human identifier BioMart table")
        mapping_path = alternatives[0]
    identifiers = pd.read_csv(mapping_path, sep="\t", dtype=str).fillna("")
    orth_by_human = {g: f for g, f in ortholog.groupby("ensembl_gene_id")}
    lookup = defaultdict(set)
    for row in identifiers.to_dict("records"):
        if row.get("uniprotswissprot"):
            lookup[row["uniprotswissprot"]].add(row["ensembl_gene_id"])
    audit, requested = [], {}
    for accession, peptides in sorted(by_protein.items()):
        genes = lookup[accession]
        # Exact accession mapping only: gene-symbol matches are not equivalent evidence.
        row = {"uniprot": accession, "humanGene": fasta[accession]["gene"], "shortPeptideCount": len(peptides)}
        if len(genes) != 1:
            audit.append({**row, "status": "no_unique_accession_to_human_gene", "humanGenes": sorted(genes)})
            continue
        human_gene = next(iter(genes))
        ortho = orth_by_human.get(human_gene)
        valid = [] if ortho is None else ortho[ortho["mmusculus_homolog_orthology_type"] == "ortholog_one2one"].to_dict("records")
        mouse_genes = {r["mmusculus_homolog_ensembl_gene"] for r in valid if r["mmusculus_homolog_ensembl_gene"]}
        if len(mouse_genes) != 1:
            audit.append({**row, "status": "no_unique_one_to_one_mouse_ortholog", "humanEnsembl": human_gene})
            continue
        mouse_gene = next(iter(mouse_genes))
        if mouse_gene not in counts.index:
            audit.append({**row, "status": "mouse_gene_absent_from_counts", "mouseEnsembl": mouse_gene})
            continue
        row.update({"humanEnsembl": human_gene, "mouseEnsembl": mouse_gene,
                    "gene": valid[0]["mmusculus_homolog_associated_gene_name"],
                    "measurable": int((t2.loc[mouse_gene, MAIN] >= 10).sum()) >= 3,
                    "status": "mapped_waiting_sequence"})
        audit.append(row)
        requested[human_gene] = None
    (cache/"join-checkpoint.json").write_text(json.dumps({"statuses": dict(Counter(r["status"] for r in audit)),
                                                       "shortPrecursors": len(by_protein),
                                                       "requests": len(requested)}, indent=2))
    print("Join checkpoint", (cache/"join-checkpoint.json").read_text(), flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        for i, (gene, data) in enumerate(pool.map(lambda g: (g, fetch_homology(cache, g)), requested)):
            requested[gene] = data
            if i % 20 == 0:
                print("Alignments", i+1, "/", len(requested), flush=True)
    cuts = predictor(cache)
    candidates = []
    for audit_row in audit:
        if audit_row["status"] != "mapped_waiting_sequence":
            continue
        accession, human_gene, mouse_gene = audit_row["uniprot"], audit_row["humanEnsembl"], audit_row["mouseEnsembl"]
        response = requested[human_gene]
        homologies = [h for d in response.get("data", []) for h in d.get("homologies", [])
                      if h["target"]["id"] == mouse_gene and h["type"] == "ortholog_one2one"]
        if len(homologies) != 1:
            audit_row["status"] = "alignment_unavailable_or_ambiguous"
            continue
        homology = homologies[0]
        human_alignment, mouse_alignment = homology["source"]["align_seq"], homology["target"]["align_seq"]
        if human_alignment.replace("-", "") != fasta[accession]["sequence"]:
            audit_row["status"] = "current_human_alignment_differs_from_pinned_precursor"
            continue
        audit_row["status"] = "sequence_verified"
        peptides = [peptide_projection(p, fasta[accession]["sequence"], human_alignment, mouse_alignment, cuts)
                    for p in by_protein[accession]]
        e = float(effects[mouse_gene])
        leave_out = [{"omitted": sample, "log2FoldChange": float(v[mouse_gene])} for sample, v in loo.items()]
        stable = all(v["log2FoldChange"]*e > 0 for v in leave_out)
        identical = sum(p["sequenceConserved"] and p["cleavageConserved"] for p in peptides)
        cleavage = sum(p["cleavageConserved"] for p in peptides)
        values = normalized.loc[mouse_gene]
        candidates.append({
            "id": mouse_gene, "gene": audit_row["gene"], "humanGene": audit_row["humanGene"],
            "humanEnsembl": human_gene, "uniprot": accession,
            "description": fasta[accession]["header"].split(" OS=")[0].split(" ", 1)[-1],
            "samples": [{"sampleId": s, "group": next(g for g, ss in GROUPS.items() if s in ss),
                         "rawCount": int(t2.loc[mouse_gene, s]), "normalizedCount": float(values[s])} for s in SAMPLES],
            "effect": {"log2FoldChange": e, "looMin": min(v["log2FoldChange"] for v in leave_out),
                       "looMax": max(v["log2FoldChange"] for v in leave_out), "stableDirection": stable,
                       "leaveOneOut": leave_out, "cpmLog2FoldChange": float(cpm_effects[mouse_gene]),
                       "normalizationDirectionAgrees": e*float(cpm_effects[mouse_gene]) > 0,
                       "contextLog2FoldChange": float(np.log2((values[GROUPS["HFD"]].mean()+1)/(values[GROUPS["LFD"]].mean()+1)))},
            "peptides": peptides,
            "evidence": {"measurable": audit_row["measurable"], "conservedPeptideCount": identical,
                         "conservedCleavagePeptideCount": cleavage,
                         "detectedN": int((t2.loc[mouse_gene, MAIN] >= 10).sum()),
                         "medianNormalized": float(values[MAIN].median()),
                         "priorityTier": "sequence-supported" if identical and audit_row["measurable"] else "audit"},
            "limitations": ["Bulk lesion RNA cannot identify secreting cell type or mature peptide abundance.",
                            "One-to-one orthology and motif conservation do not demonstrate processing, secretion or bioactivity."],
        })
    candidates.sort(key=lambda c: (not (c["evidence"]["measurable"] and c["evidence"]["conservedPeptideCount"] > 0),
                                   not c["effect"]["stableDirection"], -abs(c["effect"]["log2FoldChange"]), c["id"]))
    default = [c for c in candidates if c["evidence"]["measurable"] and c["evidence"]["conservedPeptideCount"] > 0
               and c["effect"]["stableDirection"] and abs(c["effect"]["log2FoldChange"]) >= .5]
    metadata = (cache/"GSE266899_GSM.txt").read_text()
    gsms = {re.search(r"!Sample_title = ([^\n]+)", block).group(1).strip().split("_")[0]:
            re.search(r"\^SAMPLE = (\S+)", block).group(1) for block in re.split(r"(?=\^SAMPLE)", metadata) if block.startswith("^SAMPLE")}
    summary = {"countsGenes": len(counts), "publishedPrecursors": len({p["uniprot"] for p in published}),
               "publishedPeptideRows": len(published), "shortPeptideRows": len(short), "shortPrecursors": len(by_protein),
               "mappingStatuses": dict(Counter(r["status"] for r in audit)),
               "sequenceVerifiedPrecursors": len(candidates),
               "measurablePrecursors": sum(c["evidence"]["measurable"] for c in candidates),
               "strictEligiblePrecursors": sum(c["evidence"]["measurable"] and c["evidence"]["conservedPeptideCount"] > 0 for c in candidates),
               "defaultShortlistCount": len(default), "defaultShortlistGenes": [c["gene"] for c in default],
               "reportedLengthConventionDiscrepancyRows": sum(p["actualLength"] != p["reportedLength"] for p in published)}
    data = {"schemaVersion": 1, "generatedAt": "2026-10-20", "study": {
        "accession": "GSE266899", "groups": [{"id": g, "label": {"LFD": "Low-fat diet", "HFD": "High-fat diet",
                       "SEMA": "High-fat diet + semaglutide"}[g], "n": len(ss)} for g, ss in GROUPS.items()],
        "samples": [{"id": s, "group": next(g for g, ss in GROUPS.items() if s in ss),
                     "geoAccession": gsms[s], "sizeFactor": float(sf[s])} for s in SAMPLES],
        "normalization": "Median ratios to the across-sample geometric mean; genes positive in all 11 T2 mice; no length correction.",
        "normalizationGeneCount": normalization_genes,
        "analysisType": "Exploratory descriptive precursor-RNA comparison; no p-values or differential-expression claims.",
        "unit": "Independent male mouse pancreatic neoplastic lesion, T2 (90 days after orthotopic PaNIN lesion transplantation).",
        "counts": {"total": 18, "included": 11, "excludedT1": 7},
        "mainContrast": "HFD + semaglutide (n=3 deposited) versus HFD (n=4); LFD (n=4) is context, not weight-loss matching.",
        "depositionDiscrepancy": "GEO protocol states planned n=4 per T2 group; only three SEMA columns/GSMs are deposited and analyzed.",
    }, "candidates": candidates, "summary": summary, "mappingAudit": audit,
        "selectionRules": {"lengthMin": 5, "lengthMax": 25, "minRawCount": 10, "minDetectedMainMice": 3,
            "defaultMinAbsLog2FoldChange": .5, "defaultRequireStableDirection": True,
            "defaultRequireIdenticalPeptide": True, "cleavageSpacing": 5,
            "effectFormula": "log2[(mean normalized SEMA count + 1)/(mean normalized HFD count + 1)]",
            "leaveOneOut": "Omit each of the seven HFD/SEMA mice in turn; recompute size factors from the remaining ten T2 mice. Range is sensitivity, not CI.",
            "ranking": "Eligible before audit; same-sign leave-one-out before unstable; descending absolute descriptive effect; mouse Ensembl ID tie-break.",
            "cleavageMeaning": "Both flanking sites must be called in both full human and mouse proteins by the pinned upstream motif rule (default spacing 5 aa).",
            "identityMeaning": "Exact peptide identity without alignment gaps, separately from cleavage conservation.",
            "cpmSensitivity": "Counts per million uses all deposited gene counts; +1 pseudocount in CPM units. This is a qualitative normalization sensitivity check, not a second inferential test.",
        },
        "sources": [{"title": "Svensson Lab pro-hormone-predictor", "url": "https://github.com/Svensson-Lab/pro-hormone-predictor/tree/c478c45405adce4c62baa0719891c3c89f6edcda"},
                    {"title": "Coassolo et al. Nature 2025, Supplementary Table 2", "url": "https://www.nature.com/articles/s41586-025-08683-y"},
                    {"title": "GSE266899: independent external PDAC-context data", "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266899"},
                    {"title": "Cencioni et al. 2025", "url": "https://doi.org/10.1186/s13046-024-03263-w"},
                    {"title": "Ensembl human–mouse orthology and alignments", "url": "https://rest.ensembl.org/documentation/info/homology_species_gene_id"}]}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    tsv_rows = []
    for c in candidates:
        tsv_rows.append({"mouse_ensembl": c["id"], "gene": c["gene"], "human_gene": c["humanGene"], "uniprot": c["uniprot"],
                        "measurable": c["evidence"]["measurable"], "identical_peptides_with_conserved_cleavage": c["evidence"]["conservedPeptideCount"],
                        "conserved_cleavage_peptides": c["evidence"]["conservedCleavagePeptideCount"],
                        "log2_ratio": c["effect"]["log2FoldChange"], "loo_min": c["effect"]["looMin"], "loo_max": c["effect"]["looMax"],
                        "stable_direction": c["effect"]["stableDirection"], "cpm_log2_ratio": c["effect"]["cpmLog2FoldChange"],
                        **{f"{s['sampleId']}_raw": s["rawCount"] for s in c["samples"]},
                        **{f"{s['sampleId']}_normalized": s["normalizedCount"] for s in c["samples"]}})
    pd.DataFrame(tsv_rows).to_csv(output.with_name("precursor-evidence.tsv"), sep="\t", index=False)
    pd.DataFrame(audit).to_csv(output.with_name("mapping-audit.tsv"), sep="\t", index=False)
    input_sources = json.loads((cache/"input-manifest.json").read_text())
    urls_by_file = {s["file"]: s["url"] for s in input_sources}
    if (cache/"biomart-source-qa.json").exists():
        urls_by_file.update({s["file"]: s["url"] for s in json.loads((cache/"biomart-source-qa.json").read_text())})
    manifest = {"retrievalDate": "2026-10-20", "sourceDownloads": input_sources,
        "alignmentEndpoint": "https://rest.ensembl.org/homology/id/human/{humanEnsembl}?target_species=mus_musculus;type=orthologues;content-type=application/json",
        "inputFiles": [
        {"name": p.relative_to(cache).as_posix(), "sha256": digest(p),
         **({"url": urls_by_file[p.name]} if p.name in urls_by_file else {})}
        for p in sorted(cache.rglob("*")) if p.is_file()
        and p.suffix not in {".pyc"} and not "python" in p.parts],
        "upstreamComponent": "Hash-verified FastaParser.py motif-expression AST and get_cleavage_sites only; no upstream code redistributed.",
        "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__}
    output.with_name("peptide-cancer-provenance.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.cache, args.output)
