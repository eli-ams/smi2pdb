"""Focused regression checks. Not run during implementation at user request.

From the repository root after installation: python -m unittest discover -s tests
Fixtures supply coordinates directly; no conformer generation is performed.
"""
import unittest
from types import SimpleNamespace
from contextlib import closing
from unittest.mock import patch

import numpy as np
from rdkit import Chem

from smi2pdb.packing import (
    CellList, PackingError, attempt_results, pack_molecules, random_rotation,
    sobol_points, worker_count,
)


def molecule():
    mol = Chem.MolFromSmiles("[He]")
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    mol.AddConformer(conf)
    return SimpleNamespace(type="helium", mol_obj=mol, mw_gmol=4.0026, rotate=True)


def settings(**changes):
    values = dict(
        AVOGADRO_NUMBER=6.02214076e23, target_density_kg_m3=1000.0,
        initial_density_fraction=0.30, box_shape="cubic", candidates_per_molecule=32,
        clash_scale=0.72, boundary_padding=0.0, number_of_trials=2,
        max_box_expansions=0, box_expansion_factor=1.20,
        packing_workers=1,
    )
    values.update(changes)
    return SimpleNamespace(**values)


class PackingTests(unittest.TestCase):
    def test_worker_count_is_bounded_by_attempt_budget(self):
        self.assertEqual(worker_count(8, 3), 3)
        self.assertEqual(worker_count(1, 10), 1)
        self.assertEqual(worker_count(0, 1), 1)

    def test_early_success_closes_worker_pool(self):
        with patch("smi2pdb.packing.multiprocessing.get_context") as context:
            pool = context.return_value.Pool.return_value
            pool.imap_unordered.return_value = iter([("winner", {}), ("unused", {})])
            with closing(attempt_results((), [(0, 1), (1, 2)], 2)) as results:
                self.assertEqual(next(results)[0], "winner")
            pool.terminate.assert_called_once()
            pool.join.assert_called_once()

    def test_worker_exception_closes_pool_and_propagates(self):
        with patch("smi2pdb.packing.multiprocessing.get_context") as context:
            pool = context.return_value.Pool.return_value
            pool.imap_unordered.side_effect = RuntimeError("worker error")
            with self.assertRaisesRegex(RuntimeError, "worker error"):
                with closing(attempt_results((), [(0, 1)], 2)) as results:
                    next(results)
            pool.terminate.assert_called_once()
            pool.join.assert_called_once()

    def test_spawned_workers_return_a_valid_configuration(self):
        molecules = [molecule() for _ in range(8)]
        box, report = pack_molecules(
            molecules, settings(target_density_kg_m3=100.0, packing_workers=2),
            42, lambda _: None,
        )
        self.assertTrue(report["accepted"])
        self.assertEqual(report["packing_workers_used"], 2)
        self.assertEqual(report["attempt_budget_per_box"], 2)
        self.assertIsInstance(report["winning_attempt_seed"], int)
        positions = np.array([m.mol_obj.GetConformer().GetPositions()[0] for m in molecules])
        radius = Chem.GetPeriodicTable().GetRvdw(2)
        for i in range(len(positions)):
            for j in range(i):
                delta = positions[i] - positions[j]
                delta -= np.array(box) * np.round(delta / box)
                self.assertGreaterEqual(float(np.linalg.norm(delta)), .72 * 2 * radius)

    def test_parallel_budget_exhaustion_has_no_partial_commit(self):
        # There is room for either atom individually but not both at this scale.
        molecules = [molecule(), molecule()]
        originals = [m.mol_obj.GetConformer().GetPositions().copy() for m in molecules]
        with self.assertRaises(PackingError) as caught:
            pack_molecules(molecules, settings(
                target_density_kg_m3=1400, clash_scale=1.0, packing_workers=2,
            ), 42, lambda _: None)
        self.assertFalse(caught.exception.report["accepted"])
        self.assertEqual(caught.exception.report["packing_workers_used"], 2)
        for mol, original in zip(molecules, originals):
            np.testing.assert_array_equal(mol.mol_obj.GetConformer().GetPositions(), original)

    def test_sobol_reference_prefix_and_skip(self):
        expected = [[0, 0, 0], [.5, .5, .5], [.75, .25, .25], [.25, .75, .75]]
        np.testing.assert_allclose(sobol_points(4), expected)
        np.testing.assert_array_equal(sobol_points(2, 2), sobol_points(4)[2:])

    def test_rotation_preserves_distances(self):
        rotation = random_rotation(np.random.default_rng(12))
        np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-14)
        self.assertAlmostEqual(np.linalg.det(rotation), 1.0)

    def test_periodic_clash_across_box_boundary(self):
        cells = CellList([10, 10, 10], 1.5, .72)
        cells.add(np.array([[.2, 5, 5]]), [1.5])
        self.assertTrue(cells.clashes(np.array([[9.8, 5, 5]]), [1.5]))
        self.assertFalse(cells.clashes(np.array([[5, 5, 5]]), [1.5]))

    def test_cell_checks_match_brute_force(self):
        rng = np.random.default_rng(10)
        for box in (np.array([10., 11., 12.]), np.array([3., 4., 5.])):
            positions = rng.random((20, 3)) * box
            radii = rng.uniform(.5, 1.8, 20)
            cells = CellList(box, 1.8, .72)
            cells.add(positions, radii)
            for _ in range(30):
                point = rng.random(3) * box
                radius = float(rng.uniform(.5, 1.8))
                delta = point - positions
                delta -= box * np.round(delta / box)
                expected = bool(np.any(np.sum(delta * delta, axis=1) < (.72 * (radius + radii))**2))
                self.assertEqual(cells.clashes(np.array([point]), [radius]), expected)

    def test_success_density_reproducibility_and_all_pairs(self):
        options = settings(target_density_kg_m3=100.0)
        first = [molecule() for _ in range(8)]
        second = [molecule() for _ in range(8)]
        box, report = pack_molecules(first, options, 42, lambda _: None)
        other_box, other_report = pack_molecules(second, options, 42, lambda _: None)
        self.assertEqual(report, other_report)
        self.assertEqual(box, other_box)
        self.assertTrue(report["accepted"])
        self.assertEqual(report["box_expansions"], 0)
        self.assertAlmostEqual(report["actual_initial_density_fraction"], .30, delta=.001)
        points = np.array([m.mol_obj.GetConformer().GetPositions()[0] for m in first])
        other = np.array([m.mol_obj.GetConformer().GetPositions()[0] for m in second])
        np.testing.assert_array_equal(points, other)
        radius = Chem.GetPeriodicTable().GetRvdw(2)
        for i in range(len(points)):
            for j in range(i):
                delta = points[i] - points[j]
                delta -= np.array(box) * np.round(delta / box)
                self.assertGreaterEqual(float(np.linalg.norm(delta)), .72 * 2 * radius)

    def test_failure_keeps_input_coordinates_and_reports_density(self):
        mol = molecule()
        original = mol.mol_obj.GetConformer().GetPositions().copy()
        with self.assertRaises(PackingError) as caught:
            pack_molecules([mol], settings(target_density_kg_m3=1e9), 42, lambda _: None)
        self.assertFalse(caught.exception.report["accepted"])
        self.assertEqual(caught.exception.report["box_expansions"], 0)
        np.testing.assert_array_equal(mol.mol_obj.GetConformer().GetPositions(), original)

    def test_expansion_is_explicit_and_density_is_actual(self):
        options = settings(target_density_kg_m3=8000, max_box_expansions=1,
                           box_expansion_factor=4.0)
        box, report = pack_molecules([molecule()], options, 42, lambda _: None)
        self.assertEqual(report["box_expansions"], 1)
        self.assertLess(report["actual_initial_density_fraction"], .30)
        mass_kg = 4.0026 / options.AVOGADRO_NUMBER / 1000
        self.assertAlmostEqual(report["actual_initial_density_kg_m3"], mass_kg / (np.prod(box) * 1e-30))


if __name__ == "__main__":
    unittest.main()
