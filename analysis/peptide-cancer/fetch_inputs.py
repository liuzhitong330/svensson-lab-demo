"""Download public inputs into a private cache, then report the unfiltered universe."""
from pathlib import Path
import argparse
import concurrent.futures
import gzip
import hashlib
import io
import json
import re
import urllib.request
import urllib.parse
import zipfile
from openpyxl import load_workbook
import pandas as pd
import numpy as np

COMMIT = "c478c45405adce4c62baa0719891c3c89f6edcda"
URLS = {
    "GSE266899_raw_counts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266899/suppl/GSE266899_raw_counts.txt.gz",
    "GSE266899_GSM.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266899&targ=gsm&form=text&view=full",
    "brp-supplementary-tables.zip": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-025-08683-y/MediaObjects/41586_2025_8683_MOESM3_ESM.zip",
    "secreted.fasta": f"https://raw.githubusercontent.com/Svensson-Lab/pro-hormone-predictor/{COMMIT}/secreted.fasta",
    "FastaParser.py": f"https://raw.githubusercontent.com/Svensson-Lab/pro-hormone-predictor/{COMMIT}/FastaParser.py",
}
BIOMART = {
    "human-mouse-biomart.tsv": ["ensembl_gene_id", "mmusculus_homolog_ensembl_gene",
        "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_ensembl_peptide",
        "mmusculus_homolog_orthology_type", "mmusculus_homolog_orthology_confidence",
        "mmusculus_homolog_perc_id", "mmusculus_homolog_perc_id_r1"],
    "human-uniprot-biomart.tsv": ["ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"],
}


def fetch(path, url):
    if not path.exists():
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "ResearchReanalysis/1.0"}), timeout=180) as response:
            payload = response.read()
        path.write_bytes(payload)
    return {"file": path.name, "url": url, "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def fetch_biomart(cache, filename, attributes):
    query = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
             '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">'
             '<Dataset name="hsapiens_gene_ensembl" interface="default">'
             + "".join(f'<Attribute name="{a}"/>' for a in attributes) + "</Dataset></Query>")
    url = "https://www.ensembl.org/biomart/martservice?"+urllib.parse.urlencode({"query": query})
    path = cache/filename
    if not path.exists():
        with urllib.request.urlopen(url, timeout=180) as response:
            text = response.read().decode()
        if "Query ERROR" in text or text.lstrip().startswith("<"):
            raise RuntimeError("BioMart did not return a tab-delimited export")
        path.write_text("\t".join(attributes)+"\n"+text)
    return {"file": filename, "url": url, "query": query, "attributes": attributes,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def checkpoint(cache):
    cache.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        manifest = list(pool.map(lambda kv: fetch(cache/kv[0], kv[1]), URLS.items()))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        manifest += list(pool.map(lambda kv: fetch_biomart(cache, kv[0], kv[1]), BIOMART.items()))
    counts = pd.read_csv(cache/"GSE266899_raw_counts.txt.gz", sep="\t", index_col=0)
    a = counts.to_numpy()
    assert np.isfinite(a).all() and (a >= 0).all() and np.equal(a, np.floor(a)).all()
    records = {}
    for line in (cache/"secreted.fasta").read_text().splitlines():
        if line.startswith(">"):
            accession = line.split("|")[1]
            records[accession] = {"header": line[1:], "sequence": "",
                                  "gene": re.search(r"GN=(\S+)", line).group(1) if "GN=" in line else None}
        else:
            records[accession]["sequence"] += line.strip()
    with zipfile.ZipFile(cache/"brp-supplementary-tables.zip") as archive:
        member = next(n for n in archive.namelist() if n.endswith("Supplementary Table 2.xlsx") and not n.startswith("__MACOSX"))
        wb = load_workbook(io.BytesIO(archive.read(member)), read_only=True, data_only=True)
        rows = list(wb["Supplementary Table 2"].values)
    peptides = []
    for excel_row, row in enumerate(rows[4:], 5):
        if not isinstance(row[0], str) or not isinstance(row[1], str):
            continue
        accession = row[1].split("|")[1]
        peptides.append({"sourceRow": excel_row, "sequence": row[0].strip(), "uniprot": accession,
                         "reportedLength": row[2], "actualLength": len(row[0].strip()),
                         "humanGene": records.get(accession, {}).get("gene"),
                         "matchesPinnedPrecursor": row[0].strip() in records.get(accession, {}).get("sequence", "")})
    (cache/"published-peptides.json").write_text(json.dumps(peptides, indent=2))
    (cache/"input-manifest.json").write_text(json.dumps(manifest, indent=2))
    summary = {"geneRows": len(counts), "samples": list(counts.columns), "integerNonnegativeCounts": True,
               "rawTotalPerSample": counts.sum().to_dict(), "humanFastaProteins": len(records),
               "publishedPeptideRows": len(peptides), "publishedPrecursors": len({r["uniprot"] for r in peptides}),
               "lengthMetadataMismatchRows": sum(r["reportedLength"] != r["actualLength"] for r in peptides),
               "short5to25Rows": sum(5 <= r["actualLength"] <= 25 for r in peptides),
               "peptidesNotInPinnedPrecursor": sum(not r["matchesPinnedPrecursor"] for r in peptides)}
    (cache/"checkpoint.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    return counts, records, peptides


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    checkpoint(parser.parse_args().cache)
