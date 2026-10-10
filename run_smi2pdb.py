"""Run directly from a source checkout when dependencies are already available."""
from pathlib import Path
import sys

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
    from smi2pdb.cli import main
    raise SystemExit(main())
