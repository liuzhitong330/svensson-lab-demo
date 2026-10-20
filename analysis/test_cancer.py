"""Independent standard-library spot checks of the derived cancer JSON."""
from pathlib import Path
import json
import math
import statistics
import unittest

DATA = json.loads((Path(__file__).resolve().parents[1]/"cancer-data.json").read_text())


class SourceChecks(unittest.TestCase):
    def test_counts(self):
        self.assertEqual(sum(g["n"] for g in DATA["endpoints"]["mass_change_week22"]["groups"]), 40)
        self.assertEqual(sum(g["n"] for g in DATA["endpoints"]["tumor_day16"]["groups"]), 128)
        self.assertEqual(sum(g["n"] for g in DATA["endpoints"]["immune"]["groups"]), 64)

    def test_mass_recompute(self):
        for trace in DATA["massTraces"]:
            for grams, change in zip(trace["grams"], trace["percentChange"]):
                expected = (grams-trace["grams"][0])/trace["grams"][0]*100
                self.assertAlmostEqual(change, expected, places=10)
        first = DATA["massTraces"][0]
        self.assertEqual(first["sourceRange"], "Fig. 3B!C6:I6")
        self.assertEqual(first["grams"][0], 31.9)
        self.assertEqual(first["grams"][-1], 32.6)

    def test_summaries(self):
        for endpoint in DATA["endpoints"].values():
            for group in endpoint["groups"]:
                values = group["values"]
                self.assertEqual(group["n"], len(values))
                self.assertEqual(len(values), len(group["sourceCells"]))
                self.assertAlmostEqual(group["mean"], math.fsum(values)/len(values), places=10)
                self.assertAlmostEqual(group["median"], statistics.median(values), places=10)
                variance = math.fsum((v-group["mean"])**2 for v in values)/(len(values)-1)
                self.assertAlmostEqual(group["sd"], math.sqrt(variance), places=10)

    def test_no_cross_panel_ids(self):
        self.assertNotIn("animalId", DATA)
        for endpoint in DATA["endpoints"].values():
            for group in endpoint["groups"]:
                self.assertNotIn("animalIds", group)
        self.assertEqual(len(DATA["endpoints"]["mass_change_week22"]["groups"]), 3)
        self.assertTrue(all(not g["icb"] for g in DATA["endpoints"]["mass_change_week22"]["groups"]))

    def test_reported_experimental_unit(self):
        self.assertIn("One mouse per observation", DATA["endpoints"]["tumor_day16"]["biologicalUnit"])
        self.assertIn("individual mouse", DATA["design"]["reportedExperimentalUnit"])
        self.assertIn("panel-local", DATA["endpoints"]["tumor_day16"]["biologicalUnit"])

    def test_source_spotchecks(self):
        tumor = {g["key"]: g for g in DATA["endpoints"]["tumor_day16"]["groups"]}
        immune = {g["key"]: g for g in DATA["endpoints"]["immune"]["groups"]}
        self.assertEqual(tumor["nc"]["values"][0], 415.292)
        self.assertEqual(tumor["nc"]["sourceCells"][0], "Fig. 3G!C5")
        self.assertEqual(tumor["sema"]["values"][0], 463.89)
        self.assertEqual(immune["nc_icb"]["values"], [45.9, 43.9, 55.6])


if __name__ == "__main__":
    unittest.main()
