"""Regression/invariant tests on the public derived evidence, without private inputs."""
from pathlib import Path
import argparse
import importlib.util
import json
import math
import unittest
import numpy as np
import pandas as pd

SPEC = importlib.util.spec_from_file_location("peptide_analysis", Path(__file__).with_name("analyze.py"))
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)
DATA = None


class AnalysisTests(unittest.TestCase):
    def test_median_ratio_known_scaling(self):
        frame = pd.DataFrame({"a": [10, 20, 0], "b": [20, 40, 5], "c": [40, 80, 0]})
        factors, n = analysis.ratio_factors(frame)
        self.assertEqual(n, 2)
        np.testing.assert_allclose(factors.to_numpy(), [.5, 1, 2])

    def test_zero_geometric_mean_guard(self):
        with self.assertRaises(ValueError):
            analysis.ratio_factors(pd.DataFrame({"a": [1, 0], "b": [0, 1]}))

    def test_repeated_peptide_not_guessed(self):
        p = {"sequence": "AAAAA", "sourceRow": 5, "reportedLength": 7}
        result = analysis.peptide_projection(p, "KKAAAAAKRAAAAA", "KKAAAAAKRAAAAA", "KKAAAAAKRAAAAA", lambda s: set())
        self.assertEqual(result["status"], "ambiguous_or_missing_human_location")
        self.assertFalse(result["sequenceConserved"])

    def test_indel_is_not_positionwise_identity(self):
        p = {"sequence": "ABCDE", "sourceRow": 5, "reportedLength": 7}
        human, ha, ma = "KKABCDEKR", "KKAB-CDEKR", "KKABXCDEKR"
        result = analysis.peptide_projection(p, human, ha, ma, lambda s: {(0, 2), (len(s)-2, len(s))})
        self.assertEqual(result["mouseSequence"], "ABXCDE")
        self.assertTrue(result["alignmentHasIndel"])
        self.assertFalse(result["sequenceConserved"])

    def test_mapped_substitution_is_not_identity(self):
        p = {"sequence": "ABCDE", "sourceRow": 5, "reportedLength": 7}
        result = analysis.peptide_projection(p, "KKABCDEKR", "KKABCDEKR", "KKABXDEKR",
                                            lambda s: {(0, 2), (7, 9)})
        self.assertEqual(result["identicalResidues"], 4)
        self.assertEqual(result["mismatchPositions"], [3])
        self.assertTrue(result["cleavageConserved"])
        self.assertFalse(result["sequenceConserved"])

    def test_study_independent_mouse_counts(self):
        self.assertEqual(DATA["summary"]["countsGenes"], 54309)
        self.assertEqual(DATA["study"]["normalizationGeneCount"], 13796)
        self.assertEqual({g["id"]: g["n"] for g in DATA["study"]["groups"]}, {"LFD": 4, "HFD": 4, "SEMA": 3})
        self.assertEqual(len(DATA["study"]["samples"]), 11)
        self.assertEqual(DATA["summary"]["shortPrecursors"], sum(DATA["summary"]["mappingStatuses"].values()))

    def test_all_candidate_arithmetic_and_rules(self):
        sf = {s["id"]: s["sizeFactor"] for s in DATA["study"]["samples"]}
        default = []
        for c in DATA["candidates"]:
            self.assertEqual(len(c["samples"]), 11)
            means = {}
            for group in ("LFD", "HFD", "SEMA"):
                values = [s["normalizedCount"] for s in c["samples"] if s["group"] == group]
                means[group] = sum(values)/len(values)
            observed = math.log2((means["SEMA"]+1)/(means["HFD"]+1))
            self.assertAlmostEqual(observed, c["effect"]["log2FoldChange"], places=10)
            for s in c["samples"]:
                self.assertEqual(s["rawCount"], int(s["rawCount"]))
                self.assertGreaterEqual(s["rawCount"], 0)
                self.assertAlmostEqual(s["normalizedCount"], s["rawCount"]/sf[s["sampleId"]], places=8)
            loo = [r["log2FoldChange"] for r in c["effect"]["leaveOneOut"]]
            self.assertEqual(len(loo), 7)
            self.assertEqual(c["effect"]["stableDirection"], all(v*observed > 0 for v in loo))
            self.assertEqual(c["effect"]["looMin"], min(loo))
            self.assertEqual(c["effect"]["looMax"], max(loo))
            strict = sum(p["sequenceConserved"] and p["cleavageConserved"] for p in c["peptides"])
            self.assertEqual(strict, c["evidence"]["conservedPeptideCount"])
            for p in c["peptides"]:
                if p.get("status") != "projected":
                    continue
                self.assertEqual(p["actualLength"], len(p["humanSequence"]))
                self.assertEqual(p["alignedHuman"].replace("-", ""), p["humanSequence"])
                self.assertEqual(p["alignedMouse"].replace("-", ""), p["mouseSequence"])
                self.assertEqual(p["sequenceConserved"], p["humanSequence"] == p["mouseSequence"] and not p["alignmentHasIndel"])
            if c["evidence"]["measurable"] and strict and c["effect"]["stableDirection"] and abs(observed) >= .5:
                default.append(c["gene"])
        self.assertEqual(default, DATA["summary"]["defaultShortlistGenes"])
        self.assertEqual(len(default), DATA["summary"]["defaultShortlistCount"])

    def test_brinp2_negative_audit_control(self):
        c = next(c for c in DATA["candidates"] if c["humanGene"] == "BRINP2")
        p = next(p for p in c["peptides"] if p["humanSequence"] == "THRILRRLFNLC")
        self.assertEqual(p["mouseSequence"], "MHRIVRRLFNLC")
        self.assertEqual(p["identicalResidues"], 10)
        self.assertTrue(p["cleavageConserved"])
        self.assertFalse(p["sequenceConserved"])
        self.assertFalse(c["evidence"]["measurable"])
        self.assertNotIn(c["gene"], DATA["summary"]["defaultShortlistGenes"])

    def test_no_private_paths_and_unique_rows(self):
        text = json.dumps(DATA)
        self.assertNotIn("/Users/", text)
        pairs = [(c["id"], c["uniprot"]) for c in DATA["candidates"]]
        self.assertEqual(len(pairs), len(set(pairs)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    args = parser.parse_args()
    DATA = json.loads(args.data.read_text())
    unittest.main(argv=["test_analysis.py"], verbosity=2)
