"""Configuration/entry-point checks. Not run during package preparation."""
from contextlib import redirect_stdout, redirect_stderr
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import warnings

from smi2pdb import __version__
from smi2pdb.cli import main
from smi2pdb.config import load_settings


class CliTests(unittest.TestCase):
    def test_version(self):
        output = io.StringIO()
        with redirect_stdout(output), self.assertRaises(SystemExit) as caught:
            main(["--version"])
        self.assertEqual(caught.exception.code, 0)
        self.assertIn(__version__, output.getvalue())

    def test_template_is_available_and_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "inputs.py"
            with redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--init", str(path)]), 0)
            self.assertEqual(load_settings(path).run_mode, "targeted")
            path.write_text("# my existing file\n", encoding="utf-8")
            with redirect_stderr(io.StringIO()):
                self.assertEqual(main(["--init", str(path)]), 1)
            self.assertEqual(path.read_text(), "# my existing file\n")

    def test_old_inputs_and_relative_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "old.py"
            path.write_text(
                "density = 300.0\nnumber_of_mixtures_needed = 2\nnbr_of_trials = 3\n"
                "layer_offset = 4.0\nenergy_threshold = 0.0\n"
                "molecules_dictionary = {'water': {'smiles': 'O', 'nbr_of_molecules': 2, 'rotate': True}}\n",
                encoding="utf-8",
            )
            with warnings.catch_warnings(record=True) as messages:
                warnings.simplefilter("always")
                settings = load_settings(path)
            self.assertEqual(settings.target_density_kg_m3, 300.0)
            self.assertEqual(settings.initial_density_fraction, 1.0)
            self.assertEqual(settings.mixtures_needed, 2)
            self.assertEqual(settings.number_of_trials, 3)
            self.assertEqual(settings.boundary_padding, 2.0)
            self.assertEqual(settings.molecules_dictionary["water"]["nbr_of_mols"], 2)
            self.assertEqual(Path(settings.output_dir), Path(tmp) / "output")
            self.assertEqual(len(messages), 3)

    def test_overrides_reach_runner_without_packing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "inputs.py"
            path.write_text("initial_density_fraction = 0.5\n", encoding="utf-8")
            with patch("smi2pdb.runner.main") as run:
                self.assertEqual(main(["--config", str(path), "--mode", "density_search", "--workers", "2"]), 0)
            settings = run.call_args.kwargs["settings"]
            self.assertEqual(settings.run_mode, "density_search")
            self.assertEqual(settings.packing_workers, 2)
            self.assertEqual(settings.initial_density_fraction, .5)

    def test_missing_input_returns_failure(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(io.StringIO()):
            self.assertEqual(main(["--config", str(Path(tmp) / "missing.py")]), 1)


if __name__ == "__main__":
    unittest.main()
