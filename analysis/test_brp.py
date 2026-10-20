"""Independent arithmetic and structural checks on generated brp-data.json."""
from pathlib import Path
import json
import math
import statistics
import unittest

DATA = json.loads((Path(__file__).resolve().parents[1] / "brp-data.json").read_text())
CONDITIONS = {c["id"]: c for c in DATA["conditions"]}


def choose(threshold, time_h, both_cell_lines=True):
    accepted = []
    for c in DATA["conditions"][2:]:
        assays = ["ns1", "neuro2a"] if both_cell_lines else ["neuro2a"]
        if any(c[k]["missing_count"] for k in assays):
            continue
        feeding = next(x for x in c["feeding"] if x["time_h"] == time_h)
        values = [c[k]["retained_activity"] for k in assays] + [feeding["retained_inhibition"]]
        if all(v is not None and abs(v) <= threshold for v in values):
            accepted.append(c["id"])
    return accepted


class PeptideChecks(unittest.TestCase):
    def test_source_coverage(self):
        self.assertEqual(len(CONDITIONS), 14)
        self.assertEqual(sum(c["ns1"]["n"] for c in CONDITIONS.values()), 41)
        self.assertEqual(sum(c["neuro2a"]["n"] for c in CONDITIONS.values()), 42)
        self.assertEqual(sum(len(c["feeding"]) for c in CONDITIONS.values()), 84)
        self.assertEqual(CONDITIONS["R3A"]["ns1"]["missing_count"], 1)
        self.assertEqual(CONDITIONS["R3A"]["ns1"]["n"], 2)

    def test_l8a_independent_arithmetic(self):
        # Literal source observations, not values copied from the derived JSON.
        ns1 = (statistics.mean([4.224635927, 3.65949881, 3.74437575]) - 1) / (
            statistics.mean([12.68375157, 12.5567013, 13.4247344]) - 1)
        neuro = (statistics.mean([1.458763811, 1.12878011, 2.35758298]) - 1) / (
            statistics.mean([52.8901726, 46.7765085, 85.3212395]) - 1)
        vehicle = statistics.mean([.6, .466667, .403333])
        brp = statistics.mean([0, .033333, .083333])
        mutant = statistics.mean([.466667, .466667, .343333])
        feed = (vehicle-mutant)/(vehicle-brp)
        self.assertAlmostEqual(CONDITIONS["L8A"]["ns1"]["retained_activity"], ns1, places=12)
        self.assertAlmostEqual(CONDITIONS["L8A"]["neuro2a"]["retained_activity"], neuro, places=12)
        self.assertAlmostEqual(CONDITIONS["L8A"]["feeding"][1]["retained_inhibition"], feed, places=12)

    def test_analytic_boundary_changes_decision(self):
        self.assertEqual(choose(.25, .5), ["L8A"])
        self.assertEqual(choose(.20, .5), [])
        self.assertEqual(choose(.25, 1), ["L8A"])
        self.assertEqual(choose(.25, 1.5), [])
        self.assertEqual(choose(.25, 2), [])
        self.assertEqual(choose(.30, 3), ["L8A"])
        self.assertNotIn("R3A", choose(100, .5))

    def test_no_fake_zero_or_cross_dataset_link(self):
        for c in CONDITIONS.values():
            self.assertIsNone(c["feeding"][0]["retained_inhibition"])
            self.assertIsNone(c["feeding"][0]["intake_percent_of_vehicle"])
            for k in ("ns1", "neuro2a"):
                self.assertEqual(c[k]["n"], len(c[k]["replicates"]))
                self.assertAlmostEqual(c[k]["mean"], statistics.mean(c[k]["replicates"]), places=12)
                self.assertAlmostEqual(c[k]["sd"], statistics.stdev(c[k]["replicates"]), places=12)

    def test_repository_identity(self):
        repo = DATA["repo_identity"]
        self.assertEqual(repo["commit"], "c478c45405adce4c62baa0719891c3c89f6edcda")
        self.assertEqual(repo["input_protein_count"], 2082)
        for r in repo["results"]:
            if r["gene"] == "BRINP2":
                self.assertEqual(r["cleavage_count"], 5)
                self.assertEqual(r["brp_exact_match_start_1based"], 386)
                self.assertEqual(r["brp_exact_match_end_1based"], 397)
                self.assertEqual(r["brp_flanking_sequence"], "LFKKTHRILRRLFNLCKRCH")
            else:
                self.assertEqual(r["gene"], "GCG")
                self.assertEqual(r["cleavage_count"], 6)
                self.assertIsNone(r["brp_exact_match_start_1based"])


if __name__ == "__main__":
    unittest.main()
