"""Command-line entry point; imports scientific libraries only when packing."""
import argparse
from importlib.resources import files
import multiprocessing
from pathlib import Path
import sys

from . import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate packed molecular PDB systems from SMILES.")
    parser.add_argument("--version", action="version", version=f"smi2pdb {__version__}")
    parser.add_argument("--init", metavar="PATH", help="Write an editable input template (never overwrite).")
    parser.add_argument("--config", default="inputs.py", metavar="PATH", help="Python input file (default: inputs.py).")
    parser.add_argument("--mode", choices=("targeted", "density_search"), help="Override the configured run mode.")
    parser.add_argument("--workers", type=int, metavar="N", help="Override packing workers: 0 auto, 1 serial.")
    args = parser.parse_args(argv)
    if args.workers is not None and args.workers < 0:
        parser.error("--workers must be nonnegative")
    try:
        if args.init:
            template = files("smi2pdb").joinpath("defaults.py").read_text(encoding="utf-8")
            with Path(args.init).open("x", encoding="utf-8") as stream:
                stream.write(template)
            print(f"Created {args.init}. Edit it, then run smi2pdb --config {args.init}")
            return 0
        from .config import load_settings
        settings = load_settings(args.config)
        if args.mode is not None:
            settings.run_mode = args.mode
        if args.workers is not None:
            settings.packing_workers = args.workers
        from .runner import main as run
        multiprocessing.freeze_support()
        run(settings=settings)
        return 0
    except Exception as exc:
        print(f"smi2pdb: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
