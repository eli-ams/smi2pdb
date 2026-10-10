"""Self-contained entry point for smi2pdb 2.0.0."""
import json
from pathlib import Path
import time
import shutil
from types import SimpleNamespace

from . import defaults as inputs
from .density_search import density_fractions
from .log_functions import print_to_log
from .mixture import Mixture, validate_settings
from .packing import PackingError
from .static_functions import next_job_number, write_to_pdb

from . import __version__


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


class JobRunner:
    def __init__(self, job_id, job_output_dir, settings=inputs):
        validate_settings(settings)
        self.settings = settings
        self.job_id = str(job_id)
        self.job_output_dir = Path(job_output_dir)
        self.job_output_dir.mkdir(parents=True, exist_ok=False)
        self.logfile = self.job_output_dir / "log.txt"
        self.nbr_of_mixtures_accepted = 0

    def log(self, message):
        print_to_log(self.logfile, message)

    def search_density(self, mixture, directory, index, seed):
        """Save every success and preserve the best on the first packing failure."""
        search = {
            "status": "running", "stop_reason": None, "best": None, "steps": [],
            "start_fraction": self.settings.initial_density_fraction,
            "increment": self.settings.density_search_increment,
            "maximum_fraction": self.settings.density_search_max_fraction,
            "effective_max_box_expansions": 0,
        }
        # Copy settings without mutating the user's configuration or the templates.
        options = SimpleNamespace(**{
            name: value for name, value in vars(self.settings).items()
            if not name.startswith("_")
        })
        options.max_box_expansions = 0
        options.run_mode = "targeted"
        self.log(f"Mixture {index}: density search; box expansion disabled.")
        try:
            for step, fraction in enumerate(density_fractions(
                self.settings.initial_density_fraction,
                self.settings.density_search_increment,
                self.settings.density_search_max_fraction,
            ), start=1):
                options.initial_density_fraction = fraction
                # Non-overlapping seeds for different searches and density steps.
                step_seed = seed + (step - 1) * self.settings.mixtures_needed
                step_dir = directory / "density_steps" / f"step_{step:03d}"
                step_dir.mkdir(parents=True)
                record = {
                    "step": step, "requested_fraction": fraction, "seed": step_seed,
                    "directory": str(step_dir.relative_to(directory)),
                }
                self.log(f"Mixture {index}: trying {fraction:.6g} rho "
                         f"({fraction * self.settings.target_density_kg_m3:.6g} kg/m3).")
                try:
                    report = mixture.pack(step_seed, self.log, settings=options)
                except PackingError as exc:
                    write_json(step_dir / "packing.json", exc.report)
                    record.update(status="packing_failed", error=str(exc))
                    search["steps"].append(record)
                    search["stop_reason"] = "packing_failed"
                    if search["best"] is None:
                        search["status"] = "failed"
                        raise PackingError("Density search could not pack its starting density.", exc.report) from exc
                    search["status"] = "complete"
                    self.log(f"Search stopped at {fraction:.6g} rho; highest successful requested "
                             f"fraction: {search['best']['requested_fraction']:.6g} rho.")
                    break
                write_json(step_dir / "packing.json", report)
                write_to_pdb(mixture.molecules, mixture.box_dims, step_dir / "mixture.pdb")
                record.update(
                    status="accepted",
                    actual_initial_density_fraction=report["actual_initial_density_fraction"],
                    actual_initial_density_kg_m3=report["actual_initial_density_kg_m3"],
                )
                search["steps"].append(record)
                search["best"] = dict(record)
                # Standard output paths always hold the highest saved success.
                pdb_dir = directory / "pdbs"
                pdb_dir.mkdir(exist_ok=True)
                shutil.copyfile(step_dir / "mixture.pdb", pdb_dir / f"mixture_{index}.pdb")
                shutil.copyfile(step_dir / "packing.json", directory / "packing.json")
                write_json(directory / "density_search.json", search)
            else:
                search["status"] = "complete"
                search["stop_reason"] = "reached_maximum_fraction"
                self.log(f"Mixture {index}: search reached its requested maximum fraction.")
        except Exception as exc:
            search["status"] = "failed"
            search["error"] = str(exc)
            if search["stop_reason"] is None:
                search["stop_reason"] = "error"
            raise
        finally:
            write_json(directory / "density_search.json", search)
        return {
            "index": index, "status": "accepted", "seed": search["best"]["seed"],
            "actual_initial_density_fraction": search["best"]["actual_initial_density_fraction"],
            "actual_initial_density_kg_m3": search["best"]["actual_initial_density_kg_m3"],
            "highest_successful_requested_fraction": search["best"]["requested_fraction"],
            "box_expansions": 0, "search_stop_reason": search["stop_reason"],
        }

    def run(self):
        started = time.monotonic()
        summary = {
            "version": __version__, "job_id": self.job_id, "status": "running",
            "run_mode": self.settings.run_mode,
            "requested_mixtures": self.settings.mixtures_needed,
            "accepted_mixtures": 0, "mixtures": [],
        }
        settings_names = (
            "run_mode", "density_search_increment", "density_search_max_fraction",
            "mixtures_needed", "target_density_kg_m3", "initial_density_fraction", "box_shape",
            "random_seed", "number_of_trials", "packing_workers", "candidates_per_molecule", "clash_scale",
            "boundary_padding", "max_box_expansions", "box_expansion_factor",
            "molecules_dictionary", "AVOGADRO_NUMBER",
        )
        try:
            write_json(self.job_output_dir / "settings.json", {
                name: getattr(self.settings, name) for name in settings_names
            })
            self.log("Preparing molecular conformers once for this job.")
            mixture = Mixture(self.settings)
            for index in range(1, self.settings.mixtures_needed + 1):
                seed = self.settings.random_seed + index - 1
                directory = self.job_output_dir / f"mixture_{index}"
                directory.mkdir()
                if self.settings.run_mode == "density_search":
                    result = self.search_density(mixture, directory, index, seed)
                    self.nbr_of_mixtures_accepted += 1
                    summary["accepted_mixtures"] = self.nbr_of_mixtures_accepted
                    summary["mixtures"].append(result)
                    continue
                self.log(f"Mixture {index}: packing with seed {seed}.")
                try:
                    report = mixture.pack(seed, self.log)
                except PackingError as exc:
                    write_json(directory / "packing.json", exc.report)
                    summary["mixtures"].append({"index": index, "status": "packing_failed"})
                    raise
                write_json(directory / "packing.json", report)
                pdb_dir = directory / "pdbs"
                pdb_dir.mkdir()
                write_to_pdb(mixture.molecules, mixture.box_dims, pdb_dir / f"mixture_{index}.pdb")
                self.nbr_of_mixtures_accepted += 1
                summary["accepted_mixtures"] = self.nbr_of_mixtures_accepted
                summary["mixtures"].append({
                    "index": index, "status": "accepted", "seed": seed,
                    "actual_initial_density_kg_m3": report["actual_initial_density_kg_m3"],
                    "actual_initial_density_fraction": report["actual_initial_density_fraction"],
                    "box_expansions": report["box_expansions"],
                })
                self.log(f"Mixture {index} accepted at {report['actual_initial_density_kg_m3']:.6g} kg/m3 "
                         f"({report['actual_initial_density_fraction']:.3%} of target); "
                         f"{report['box_expansions']} box expansions.")
            summary["status"] = "complete"
            self.log(f"Completed all {self.nbr_of_mixtures_accepted} requested mixtures.")
        except Exception as exc:
            summary["status"] = "failed"
            summary["error"] = str(exc)
            self.log(f"Job failed after {self.nbr_of_mixtures_accepted} accepted mixtures: {exc}")
            raise
        finally:
            summary["elapsed_seconds"] = time.monotonic() - started
            write_json(self.job_output_dir / "job_summary.json", summary)


def main(job_id=None, job_output_dir=None, *, settings=inputs):
    validate_settings(settings)
    if job_output_dir is None:
        root = Path(settings.output_dir)
        if not root.is_absolute():
            root = Path.cwd() / root
        job_id = str(next_job_number(root)) if job_id is None else str(job_id)
        job_output_dir = root / job_id
    else:
        job_id = str(job_id) if job_id is not None else Path(job_output_dir).name
    JobRunner(job_id, job_output_dir, settings).run()


if __name__ == "__main__":
    main()
