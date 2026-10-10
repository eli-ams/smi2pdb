"""Load user-owned Python configuration without importing a global 'inputs'."""
from copy import deepcopy
from pathlib import Path
import runpy
from types import SimpleNamespace
import warnings

from . import defaults


def load_settings(path):
    """Execute a trusted local inputs file and merge its values over defaults."""
    path = Path(path).resolve()
    supplied = runpy.run_path(str(path))
    values = {name: deepcopy(value) for name, value in vars(defaults).items()
              if not name.startswith("_")}
    for old, new in (("number_of_mixtures_needed", "mixtures_needed"),
                     ("nbr_of_trials", "number_of_trials")):
        if old in supplied and new not in supplied:
            supplied[new] = supplied[old]
    if "target_density_kg_m3" not in supplied:
        legacy_density = supplied.get("density", supplied.get("initial_density"))
        if legacy_density is not None:
            supplied["target_density_kg_m3"] = legacy_density
            supplied["initial_density_fraction"] = 1.0
            warnings.warn("Legacy density is treated as an absolute packing density in kg/m3.", stacklevel=2)
    if "layer_offset" in supplied and "boundary_padding" not in supplied:
        supplied["boundary_padding"] = supplied["layer_offset"] / 2
        warnings.warn("layer_offset maps to twice the boundary padding; molecular extent margins also apply.", stacklevel=2)
    if "energy_threshold" in supplied or "cutoff_distance" in supplied:
        warnings.warn("Legacy energy_threshold/cutoff_distance are ignored: acceptance now uses vdW distances.", stacklevel=2)
    values.update({name: supplied[name] for name in values if name in supplied})
    values = deepcopy(values)
    if not isinstance(values["molecules_dictionary"], dict):
        raise ValueError("molecules_dictionary must be a dictionary.")
    for info in values["molecules_dictionary"].values():
        if isinstance(info, dict) and "nbr_of_molecules" in info and "nbr_of_mols" not in info:
            info["nbr_of_mols"] = info.pop("nbr_of_molecules")
    output = Path(values["output_dir"]).expanduser()
    values["output_dir"] = str(output if output.is_absolute() else path.parent / output)
    return SimpleNamespace(**values)
