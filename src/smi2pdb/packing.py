"""Incremental Sobol packing, adapted from the scymol reference implementation.

All runtime dependencies are local modules, NumPy, and RDKit. Intermolecular
contacts are checked against scaled vdW radii under periodic boundaries.
"""

import math
import multiprocessing
import os
from contextlib import closing
from itertools import product

import numpy as np
from rdkit import Chem
from . import __version__


def sobol_points(count, skip=0):
    """Deterministic three-dimensional Sobol sequence (32-bit direction numbers)."""
    if count < 0 or skip < 0 or skip + count > 2**32:
        raise ValueError("Sobol indices must be within the unsigned 32-bit range.")
    directions = np.zeros((3, 32), dtype=np.uint32)
    directions[0] = [np.uint32(1 << (31 - i)) for i in range(32)]
    for dimension, (degree, coefficient, initial) in enumerate(
        ((1, 0, (1,)), (2, 1, (1, 3))), start=1
    ):
        for i, value in enumerate(initial):
            directions[dimension, i] = np.uint32(value << (31 - i))
        for i in range(degree, 32):
            value = directions[dimension, i - degree]
            value ^= value >> np.uint32(degree)
            for offset in range(1, degree):
                if (coefficient >> (degree - 1 - offset)) & 1:
                    value ^= directions[dimension, i - offset]
            directions[dimension, i] = value
    points = np.empty((count, 3), dtype=float)
    for output_index, sequence_index in enumerate(range(skip, skip + count)):
        gray = sequence_index ^ (sequence_index >> 1)
        value = np.zeros(3, dtype=np.uint32)
        bit = 0
        while gray:
            if gray & 1:
                value ^= directions[:, bit]
            gray >>= 1
            bit += 1
        points[output_index] = value.astype(float) / float(2**32)
    return points


def random_rotation(rng):
    """Uniform rotation sampled through a unit quaternion."""
    u1, u2, u3 = rng.random(3)
    x = math.sqrt(1 - u1) * math.sin(2 * math.pi * u2)
    y = math.sqrt(1 - u1) * math.cos(2 * math.pi * u2)
    z = math.sqrt(u1) * math.sin(2 * math.pi * u3)
    w = math.sqrt(u1) * math.cos(2 * math.pi * u3)
    return np.array([
        [1 - 2*(y*y + z*z), 2*(x*y - z*w), 2*(x*z + y*w)],
        [2*(x*y + z*w), 1 - 2*(x*x + z*z), 2*(y*z - x*w)],
        [2*(x*z - y*w), 2*(y*z + x*w), 1 - 2*(x*x + y*y)],
    ])


class CellList:
    """Periodic spatial bins containing only previously accepted molecules."""

    def __init__(self, box, maximum_radius, scale):
        self.box = np.asarray(box, dtype=float)
        self.scale = scale
        # Every bin is at least as wide as the largest possible clash distance.
        cell_size = max(2.5, 2 * maximum_radius * scale)
        self.counts = np.maximum(1, np.floor(self.box / cell_size).astype(int))
        self.widths = self.box / self.counts
        self.bins = {}

    def key(self, position):
        return tuple((np.floor(position / self.widths).astype(int) % self.counts).tolist())

    def clashes(self, coordinates, radii):
        for position, radius in zip(coordinates, radii):
            cell = self.key(position)
            neighbours = {
                tuple((cell[axis] + offset[axis]) % int(self.counts[axis])
                      for axis in range(3))
                for offset in product((-1, 0, 1), repeat=3)
            }
            for neighbour in neighbours:
                for other_position, other_radius in self.bins.get(neighbour, ()):
                    delta = position - other_position
                    delta -= self.box * np.round(delta / self.box)
                    limit = self.scale * (radius + other_radius)
                    if np.dot(delta, delta) < limit * limit:
                        return True
        return False

    def add(self, coordinates, radii):
        for position, radius in zip(coordinates, radii):
            self.bins.setdefault(self.key(position), []).append((position.copy(), float(radius)))


class PackingError(RuntimeError):
    """No full configuration found within the configured search budget."""

    def __init__(self, message, report):
        super().__init__(message)
        self.report = report


def worker_count(requested, trials):
    """Zero selects available CPUs minus one, bounded by the attempt budget."""
    available = getattr(os, "process_cpu_count", os.cpu_count)() or 1
    count = max(1, available - 1) if requested == 0 else requested
    return min(count, trials)


def _packing_attempt(payload, task):
    """One independent attempt; only arrays/scalars cross the process boundary."""
    insertion_order, box, maximum_radius, scale, padding, candidate_count = payload
    attempt, attempt_seed = task
    rng = np.random.default_rng(attempt_seed)
    skip = int(rng.integers(1, 100_001))
    candidates = (sobol_points(candidate_count, skip) + rng.random(3)) % 1.0
    cells = CellList(box, maximum_radius, scale)
    placed = {}
    candidate_index = 0
    rejected = 0
    reason = "candidate budget exhausted"
    for index, name, rotate, base, radii, extent in insertion_order:
        margin = extent + padding + 0.001
        usable = box - 2 * margin
        if np.any(usable <= 0):
            reason = f"Box too small for {name}'s bounding extent"
            break
        while candidate_index < candidate_count:
            centre = margin + candidates[candidate_index] * usable
            candidate_index += 1
            rotation = random_rotation(rng) if rotate else np.eye(3)
            coordinates = np.round(base @ rotation.T + centre, 3)
            if cells.clashes(coordinates, radii):
                rejected += 1
                continue
            cells.add(coordinates, radii)
            placed[index] = coordinates
            break
        if index not in placed:
            break
    complete = len(placed) == len(insertion_order)
    record = {
        "attempt": attempt + 1, "attempt_seed": attempt_seed, "sobol_skip": skip,
        "box_dimensions_angstrom": box.tolist(), "placed_molecules": len(placed),
        "candidates_tested": candidate_index, "rejected_candidates": rejected,
        "accepted": complete, "reason": "complete" if complete else reason,
    }
    return (placed if complete else None), record


_worker_payload = None


def _initialize_worker(payload):
    global _worker_payload
    _worker_payload = payload


def _worker_attempt(task):
    return _packing_attempt(_worker_payload, task)


def attempt_results(payload, tasks, workers):
    """Yield completed attempts; closing terminates and joins owned workers.

The parent alone writes output. Explicit spawn works on Windows and avoids
forking a process that has already initialized RDKit/native libraries.
"""
    if workers == 1:
        for task in tasks:
            yield _packing_attempt(payload, task)
        return
    pool = multiprocessing.get_context("spawn").Pool(
        processes=workers, initializer=_initialize_worker, initargs=(payload,)
    )
    try:
        yield from pool.imap_unordered(_worker_attempt, tasks, chunksize=1)
    finally:
        # Also runs on first success, an exception, or keyboard interruption.
        pool.terminate()
        pool.join()


def pack_molecules(molecules, settings, seed, progress=print):
    """Pack rigid conformers; commit coordinates only after a complete success.

Coordinates are rounded before checking, so the exported PDB coordinates meet
the same distance criterion. Extent margins also keep molecules clear of their
own periodic images. Existing placements are preserved within each attempt.
"""
    table = Chem.GetPeriodicTable()
    templates = []
    for index, molecule in enumerate(molecules):
        coordinates = np.array(molecule.mol_obj.GetConformer().GetPositions(), copy=True)
        coordinates -= coordinates.mean(axis=0)
        radii = np.array([table.GetRvdw(a.GetAtomicNum()) for a in molecule.mol_obj.GetAtoms()])
        if not np.all(np.isfinite(coordinates)) or not np.all(np.isfinite(radii)) or np.any(radii <= 0):
            raise ValueError(f"Invalid geometry or vdW radii for {molecule.type}.")
        extent = float(np.max(np.linalg.norm(coordinates, axis=1) + radii))
        templates.append((index, molecule.type, molecule.rotate, coordinates, radii, extent))
    insertion_order = sorted(templates, key=lambda item: item[5], reverse=True)
    total_mass = sum(m.mw_gmol for m in molecules)
    target_volume = total_mass / settings.AVOGADRO_NUMBER / settings.target_density_kg_m3 * 1e27
    final_length = target_volume ** (1 / 3)
    fraction = settings.initial_density_fraction
    if settings.box_shape == "cubic":
        box = np.full(3, (target_volume / fraction) ** (1 / 3))
    else:
        box = np.array([final_length, final_length, final_length / fraction])
    # PDB CRYST1 stores three decimal places. Round upwards before packing.
    box = np.ceil(box * 1000) / 1000
    maximum_radius = max(float(np.max(item[4])) for item in templates)
    candidate_count = max(512, len(molecules) * settings.candidates_per_molecule)
    requested_workers = getattr(settings, "packing_workers", 1)
    workers = worker_count(requested_workers, settings.number_of_trials)
    history = []
    accepted = None
    winning_attempt_seed = None
    for expansion in range(settings.max_box_expansions + 1):
        tasks = [
            (attempt, int(np.random.SeedSequence([seed, expansion, attempt]).generate_state(
                1, dtype=np.uint64
            )[0]))
            for attempt in range(settings.number_of_trials)
        ]
        payload = (insertion_order, box, maximum_radius, settings.clash_scale,
                   settings.boundary_padding, candidate_count)
        progress(f"Packing with {workers} worker(s), at most {settings.number_of_trials} attempts: "
                 f"box {box.tolist()} angstrom; density "
                 f"{settings.target_density_kg_m3 * target_volume / np.prod(box):.6g} kg/m3")
        with closing(attempt_results(payload, tasks, workers)) as results:
            for placed, record in results:
                record["expansion"] = expansion
                history.append(record)
                progress(f"Attempt {record['attempt']}: placed {record['placed_molecules']}/"
                         f"{len(molecules)}; rejected {record['rejected_candidates']} candidates.")
                if placed is not None:
                    accepted = placed
                    winning_attempt_seed = record["attempt_seed"]
                    break
                if record["reason"].startswith("Box too small"):
                    break
        if accepted is not None:
            break
        if expansion < settings.max_box_expansions:
            box = np.ceil(box * settings.box_expansion_factor * 1000) / 1000
            progress("Expanding the box; the actual density will be lower than requested.")
    actual_fraction = float(target_volume / np.prod(box))
    report = {
        "version": __version__, "method": "sobol-vdw-rejection", "accepted": accepted is not None,
        "random_seed": seed, "molecule_count": len(molecules), "box_shape": settings.box_shape,
        "packing_workers_requested": requested_workers, "packing_workers_used": workers,
        "attempt_budget_per_box": settings.number_of_trials,
        "winning_attempt_seed": winning_attempt_seed,
        "attempt_statistics_scope": "results received before stopping; cancelled work excluded",
        "target_density_kg_m3": settings.target_density_kg_m3,
        "requested_initial_density_fraction": fraction,
        "actual_initial_density_fraction": actual_fraction,
        "actual_initial_density_kg_m3": settings.target_density_kg_m3 * actual_fraction,
        "box_dimensions_angstrom": box.tolist(), "box_expansions": expansion,
        "clash_scale": settings.clash_scale, "boundary_padding_angstrom": settings.boundary_padding,
        "candidate_budget_per_attempt": candidate_count,
        "rejected_candidates": sum(h["rejected_candidates"] for h in history), "attempts": history,
    }
    if accepted is None:
        raise PackingError("Could not place every molecule within the packing budget.", report)
    for index, molecule in enumerate(molecules):
        conformer = molecule.mol_obj.GetConformer()
        for atom_index, position in enumerate(accepted[index]):
            conformer.SetAtomPosition(atom_index, tuple(float(v) for v in position))
    return box.tolist(), report
