"""Reanalyse published peptide-control data without modifying source workbooks.

Usage: python analysis/brp_analysis.py --input-dir /path/to/inputs --output brp-data.json
Dependencies: openpyxl, regex.
Fetch the pinned repository and source workbook described in source_manifest_brp.json.
No upstream repository code is redistributed by this analysis.
"""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import math
import statistics
import sys
from types import SimpleNamespace

ARGP = argparse.ArgumentParser()
ARGP.add_argument("--input-dir", type=Path, required=True)
ARGP.add_argument("--output", type=Path, default=Path("brp-data.json"))
ARGS = ARGP.parse_args()
ROOT = ARGS.input_dir.resolve()
sys.path.insert(0, str(ROOT / "python-deps"))
import regex
from openpyxl import load_workbook

REPO = ROOT / "pro-hormone-predictor"
REPO_COMMIT = "c478c45405adce4c62baa0719891c3c89f6edcda"
XLSX = ROOT / "brp-source-edfig10.xlsx"
BRP = "THRILRRLFNLC"
EXPECTED_SHA256 = {
    "brp-source-edfig10.xlsx": "fab39e26f8650ea0216e233cb9e250b7642dcc3ef60a9421d1f997f859fe41d0",
    "secreted.fasta": "35cad8ba94f0a201f706826a789adabfc744023132f265d368b27988875be6d1",
    "FastaParser.py": "1ad36e296cf66e7d488fed5cdd765cb1c8ddb07bb27d4c8abff66527e4cf5436",
}


def verify_source(path):
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != EXPECTED_SHA256[path.name]:
        raise ValueError(f"Source checksum mismatch: {path.name}. Inspect any source revision before rerunning.")


def stats(values):
    valid = [float(v) for v in values if isinstance(v, (int, float))]
    return {
        "replicates": valid,
        "n": len(valid),
        "missing_count": len(values) - len(valid),
        "mean": statistics.mean(valid),
        "sd": statistics.stdev(valid) if len(valid) > 1 else None,
    }


def fasta_records(path):
    records = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append({"header": line[1:], "sequence": ""})
        else:
            records[-1]["sequence"] += line.strip()
    return records


def run_repo_component():
    """Execute the inspected upstream cleavage-pattern expression and method.

    Do not run its network-calling tissue_specificity or broken FASTA writer.
    An adapter supplies sequence_as_string(), avoiding irrelevant fastaparser.
    """
    verify_source(REPO / "FastaParser.py")
    verify_source(REPO / "secreted.fasta")
    source = (REPO / "FastaParser.py").read_text()
    tree = ast.parse(source)
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef))
    init = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "__init__")
    assignment = next(
        n for n in init.body
        if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Attribute) and t.attr == "prohormone_regex" for t in n.targets)
    )
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "get_cleavage_sites")
    namespace = {"regex": regex}
    exec(compile(ast.Module(body=[method], type_ignores=[]), "upstream-get_cleavage_sites", "exec"), namespace)
    constants = SimpleNamespace(SIGNAL_PEPTIDE_LENGTH=20)
    records = fasta_records(REPO / "secreted.fasta")
    results = []
    for gene in ("BRINP2", "GCG"):
        record = next(r for r in records if f"GN={gene} " in r["header"])
        for spacing in (5, 3):
            instance = SimpleNamespace(cleavage_spacing=spacing)
            instance.prohormone_regex = eval(
                compile(ast.Expression(assignment.value), "upstream-regex-expression", "eval"),
                {"constants": constants, "self": instance},
            )
            sequence = record["sequence"]
            adapter = SimpleNamespace(sequence_as_string=lambda: sequence)
            matches = namespace["get_cleavage_sites"](instance, adapter)
            cuts = [{"start_1based": m.start()+1, "end_1based": m.end(), "motif": m.group()} for m in matches]
            brp_start = sequence.find(BRP)
            results.append({
                "gene": gene, "uniprot": record["header"].split("|")[1],
                "length_aa": len(sequence), "minimum_spacing_aa": spacing,
                "cleavage_sites": cuts, "cleavage_count": len(cuts),
                "passes_at_least_four_sites": len(cuts) >= 4,
                "brp_exact_match_start_1based": brp_start+1 if brp_start >= 0 else None,
                "brp_exact_match_end_1based": brp_start+len(BRP) if brp_start >= 0 else None,
                "brp_flanking_sequence": sequence[max(0,brp_start-4):brp_start+len(BRP)+4] if brp_start>=0 else None,
            })
    return {
        "repository": "https://github.com/Svensson-Lab/pro-hormone-predictor",
        "commit": REPO_COMMIT,
        "license": "No LICENSE file at pinned commit; upstream code not redistributed.",
        "fasta_sha256": hashlib.sha256((REPO/"secreted.fasta").read_bytes()).hexdigest(),
        "parser_sha256": hashlib.sha256((REPO/"FastaParser.py").read_bytes()).hexdigest(),
        "execution": "Executed upstream get_cleavage_sites method and upstream regex-expression AST over checked-in FASTA, via a sequence adapter. Did not run full main.py, tissue_specificity, or FASTA writer.",
        "limitations": [
            "A motif-based precursor prediction is not evidence of secretion, processing, amidation, receptor engagement, or biological activity.",
            "Default repository spacing is 5 aa; the paper describes a 3-aa spacing screen. Both are inspected, not silently harmonized.",
            "The inspected upstream split_fasta_sequence writer appends trailing fragments inside its loop and has a for/else append; it was not relied upon.",
        ],
        "input_protein_count": len(records),
        "results": results,
    }


def main():
    verify_source(XLSX)
    wb = load_workbook(XLSX, read_only=True, data_only=True)
    ns1 = list(wb["Extended Data Fig 10a"].values)
    n2a = list(wb["Extended Data Fig 10b"].values)
    feeding = list(wb["Extended Data Fig 10c"].values)
    assert ns1[4] == n2a[4]
    conditions = []
    for j, label in enumerate(ns1[4]):
        seq = BRP if label == "BRP" else None if label == "Vehicle" else label
        c = {
            "id": "vehicle" if label == "Vehicle" else "brp" if label == "BRP" else f"{BRP[j-2]}{j-1}A",
            "label": label, "sequence": seq,
            "mutation_position": j-1 if j>=2 else None,
            "ns1": stats([ns1[r][j] for r in (5,6,7)]),
            "neuro2a": stats([n2a[r][j] for r in (5,6,7)]),
            "feeding": [],
        }
        if seq and label not in ("BRP", "Vehicle"):
            changed = [i+1 for i,(a,b) in enumerate(zip(BRP,seq)) if a != b]
            assert changed == [j-1] and seq[j-2] == "A"
        col = 1 + j*3
        annotated = feeding[4][col]
        assert regex.sub(r"\d+", "", annotated) == label
        for row in feeding[5:]:
            c["feeding"].append({
                "time_h": row[0],
                **stats(list(row[col:col+3])),
            })
        conditions.append(c)
    vehicle, brp = conditions[:2]
    for c in conditions:
        for assay in ("ns1", "neuro2a"):
            denom = brp[assay]["mean"] - vehicle[assay]["mean"]
            c[assay]["retained_activity"] = (c[assay]["mean"]-vehicle[assay]["mean"])/denom
            c[assay]["retained_replicates"] = [(v-vehicle[assay]["mean"])/denom for v in c[assay]["replicates"]]
        for i, d in enumerate(c["feeding"]):
            ref, parent = vehicle["feeding"][i]["mean"], brp["feeding"][i]["mean"]
            d["retained_inhibition"] = (ref-d["mean"])/(ref-parent) if ref != parent else None
            d["intake_percent_of_vehicle"] = 100*d["mean"]/ref if ref else None
            d["retained_experiment_values"] = [(ref-v)/(ref-parent) for v in d["replicates"]] if ref != parent else None
    provenance = {
        "publication": "Coassolo et al. Nature 2025, 641:192–201",
        "doi": "10.1038/s41586-025-08683-y",
        "source_url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-025-08683-y/MediaObjects/41586_2025_8683_MOESM19_ESM.xlsx",
        "sha256": hashlib.sha256(XLSX.read_bytes()).hexdigest(),
        "sheets": ["Extended Data Fig 10a", "Extended Data Fig 10b", "Extended Data Fig 10c"],
        "cell_rows": "Headers row 5, biological replicates rows 6:8; columns A:N",
        "feeding_rows": "Headers row 5, time and experiment means rows 6:11; columns A:AQ",
        "cell_unit": "Fos expression relative to vehicle, 1 h treatment, 100 micrograms/ml; 3 biological replicates per condition except NS-1 R3A has one technical error.",
        "feeding_unit": "Cumulative food intake in grams, 0–3 h after 5 mg/kg IP in lean male C57BL/6J mice; 3 independent experiment means, 3 mice/group per experiment; individual-animal values not supplied.",
        "normalization": "Cell retained activity = (condition mean − vehicle mean)/(BRP mean − vehicle mean). Feeding retained inhibition = (vehicle mean − condition mean)/(vehicle mean − BRP mean), separately for each time point. Reference means are not paired with individual replicate IDs.",
        "uncertainty": "Show all reported biological replicate values (cells) or independent experiment means (feeding) and sample SD. No significance, equivalence, or biological replication count is inferred from pseudo-paired columns.",
        "reference_normalization": "Vehicle Fos values fixed at 1 are source-normalized values, not evidence of zero raw assay noise. Normalized replicate spreads condition on the group reference means; they are not confidence intervals.",
        "summary_sheet_not_used": "Recomputed from Extended Data 10a–c, not 10d. The 10d feeding percentages use different aggregation from the ratio of group means used here. Its BRP Neuro2a summary is 58.3293069, whereas the mean of its three explicit 10b observations is 61.6626402. The cause of this summary difference is not established.",
        "reuse": "Paper/source data are publicly accessible but article has publisher-exclusive rights. Publish only minimal derived numerical facts and original transformation code, with source citation; do not republish workbooks, figures, or upstream code.",
        "limitation": "BRP is a non-incretin peptide. These are not semaglutide-treated tumor data, and this control audit cannot establish GLP-1 receptor engagement or cancer effects.",
    }
    output = {
        "provenance": provenance,
        "conditions": conditions,
        "repo_identity": run_repo_component(),
    }
    ARGS.output.write_text(json.dumps(output, indent=2) + "\n")
    for c in conditions[2:]:
        values = [c["ns1"]["retained_activity"], c["neuro2a"]["retained_activity"]]
        feed = next(x for x in c["feeding"] if x["time_h"] == .5)
        print(c["id"], "cell retained", [round(v,4) for v in values],
              "feeding 30min retained", round(feed["retained_inhibition"],4))
    print(json.dumps(output["repo_identity"], indent=2))


if __name__ == "__main__":
    main()

