"""Search-control regression tests; no conformer generation required.

Not executed during implementation. Run with: python -m unittest discover -s tests
"""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from smi2pdb import defaults as inputs
from smi2pdb.density_search import density_fractions
from smi2pdb.packing import PackingError
from smi2pdb.runner import JobRunner


def configuration():
    values = {k: v for k, v in vars(inputs).items() if not k.startswith("_")}
    values.update(run_mode="density_search", initial_density_fraction=.3,
                  density_search_increment=.2, density_search_max_fraction=.7,
                  max_box_expansions=5, mixtures_needed=1)
    return SimpleNamespace(**values)


class FakeMixture:
    def __init__(self, fail_at=None):
        self.fail_at = fail_at
        self.calls = []
        self.molecules = []
        self.box_dims = [10., 10., 10.]

    def pack(self, seed, progress, *, settings):
        fraction = settings.initial_density_fraction
        self.calls.append((fraction, settings.max_box_expansions, seed))
        report = {"accepted": True, "requested_initial_density_fraction": fraction,
                  "actual_initial_density_fraction": fraction - .00001,
                  "actual_initial_density_kg_m3": (fraction - .00001) * 1000,
                  "box_expansions": 0}
        if self.fail_at is not None and fraction >= self.fail_at:
            report["accepted"] = False
            raise PackingError("Budget exhausted", report)
        return report


class DensitySearchTests(unittest.TestCase):
    def test_schedule_includes_upper_endpoint(self):
        self.assertEqual(list(density_fractions(.3, .2, 1)), [.3, .5, .7, .9, 1.])
        self.assertEqual(list(density_fractions(1., .1, 1.)), [1.])

    def test_schedule_rejects_invalid_ranges(self):
        for args in ((0, .1, 1), (.3, 0, 1), (.3, .1, 1.1), (.8, .1, .7), (.3, float("nan"), 1)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                list(density_fractions(*args))

    def run_search(self, root, fake):
        settings = configuration()
        runner = JobRunner("1", root / "job", settings)
        directory = runner.job_output_dir / "mixture_1"
        directory.mkdir()
        def export(molecules, box, output):
            Path(output).write_text(str(fake.calls[-1][0]), encoding="ascii")
        with patch("smi2pdb.runner.write_to_pdb", side_effect=export):
            result = runner.search_density(fake, directory, 1, 42)
        return result, directory, settings

    def test_failure_preserves_best_and_disables_expansion(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake = FakeMixture(fail_at=.7)
            result, directory, settings = self.run_search(Path(tmp), fake)
            self.assertEqual(result["highest_successful_requested_fraction"], .5)
            self.assertEqual(result["search_stop_reason"], "packing_failed")
            self.assertEqual(fake.calls, [(.3, 0, 42), (.5, 0, 43), (.7, 0, 44)])
            self.assertEqual(settings.max_box_expansions, 5)
            self.assertEqual(settings.initial_density_fraction, .3)
            self.assertEqual((directory / "pdbs/mixture_1.pdb").read_text(), "0.5")
            report = json.loads((directory / "packing.json").read_text())
            self.assertEqual(report["requested_initial_density_fraction"], .5)
            self.assertFalse((directory / "density_steps/step_003/mixture.pdb").exists())
            search = json.loads((directory / "density_search.json").read_text())
            self.assertEqual(search["status"], "complete")
            self.assertEqual(len(search["steps"]), 3)

    def test_success_stops_at_maximum(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, directory, _ = self.run_search(Path(tmp), FakeMixture())
            self.assertEqual(result["highest_successful_requested_fraction"], .7)
            self.assertEqual(result["search_stop_reason"], "reached_maximum_fraction")
            self.assertEqual((directory / "pdbs/mixture_1.pdb").read_text(), "0.7")

    def test_starting_failure_produces_no_best_pdb(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(PackingError):
                self.run_search(Path(tmp), FakeMixture(fail_at=.3))
            directory = Path(tmp) / "job/mixture_1"
            search = json.loads((directory / "density_search.json").read_text())
            self.assertEqual(search["status"], "failed")
            self.assertIsNone(search["best"])
            self.assertFalse((directory / "pdbs/mixture_1.pdb").exists())


if __name__ == "__main__":
    unittest.main()
