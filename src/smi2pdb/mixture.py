"""Mixture configuration validation and molecular template preparation."""
import math

from . import defaults as inputs
from .molecule import Molecule
from .packing import pack_molecules


def validate_settings(settings=inputs):
    if settings.run_mode not in ("targeted", "density_search"):
        raise ValueError("run_mode must be 'targeted' or 'density_search'.")
    if settings.run_mode == "density_search":
        increment = settings.density_search_increment
        maximum = settings.density_search_max_fraction
        if not math.isfinite(increment) or not 0 < increment <= 1:
            raise ValueError("density_search_increment must be in (0, 1].")
        if not math.isfinite(maximum) or not settings.initial_density_fraction <= maximum <= 1:
            raise ValueError("Search maximum must be between initial_density_fraction and 1.0.")
    for name in ("mixtures_needed", "number_of_trials", "candidates_per_molecule"):
        value = getattr(settings, name)
        if type(value) is not int or value < 1:
            raise ValueError(f"{name} must be a positive integer.")
    for name in ("random_seed", "max_box_expansions", "packing_workers"):
        value = getattr(settings, name)
        if type(value) is not int or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer.")
    for name in ("target_density_kg_m3", "AVOGADRO_NUMBER"):
        value = getattr(settings, name)
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive.")
    for name in ("initial_density_fraction", "clash_scale"):
        value = getattr(settings, name)
        if not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError(f"{name} must be in (0, 1].")
    if not math.isfinite(settings.boundary_padding) or settings.boundary_padding < 0:
        raise ValueError("boundary_padding must be finite and nonnegative.")
    if not math.isfinite(settings.box_expansion_factor) or settings.box_expansion_factor <= 1:
        raise ValueError("box_expansion_factor must be finite and greater than 1.")
    if settings.box_shape not in ("cubic", "elongated_z"):
        raise ValueError("box_shape must be 'cubic' or 'elongated_z'.")
    if not isinstance(settings.molecules_dictionary, dict) or not settings.molecules_dictionary:
        raise ValueError("Define at least one molecular species.")
    for name, info in settings.molecules_dictionary.items():
        if not isinstance(info, dict):
            raise ValueError(f"Configuration for {name!r} must be a dictionary.")
        count = info.get("nbr_of_mols")
        if type(count) is not int or count < 1:
            raise ValueError(f"Molecule count for {name!r} must be a positive integer.")
        if not isinstance(info.get("smiles"), str) or not info["smiles"].strip():
            raise ValueError(f"Provide a SMILES string for {name!r}.")
        if type(info.get("rotate")) is not bool:
            raise ValueError(f"rotate for {name!r} must be True or False.")


class Mixture:
    def __init__(self, settings=inputs):
        validate_settings(settings)
        self.settings = settings
        self.molecules = []
        for index, (name, info) in enumerate(settings.molecules_dictionary.items()):
            seed = (settings.random_seed + index) % (2**31 - 1)
            parent = Molecule(name, info, seed)
            self.molecules.extend(parent.copy() for _ in range(info["nbr_of_mols"]))
        self.nbr_of_mols = len(self.molecules)
        self._templates = [molecule.copy() for molecule in self.molecules]
        self.box_dims = None
        self.packing = None
        self.accepted = False

    def pack(self, seed, progress=print, *, settings=None):
        self.accepted = False
        self.box_dims = None
        self.packing = None
        # Each mixture starts from the same unrounded, unplaced conformers.
        self.molecules = [molecule.copy() for molecule in self._templates]
        self.box_dims, self.packing = pack_molecules(
            self.molecules, self.settings if settings is None else settings, seed, progress
        )
        self.accepted = True
        return self.packing
