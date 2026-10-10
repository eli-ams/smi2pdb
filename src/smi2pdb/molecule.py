"""Seeded RDKit conformers used by the local packing engine."""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


class Molecule:
    def __init__(self, name, info, seed=42):
        self.type = name
        self.smiles = info["smiles"]
        self.rotate = info["rotate"]
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            raise ValueError(f"Invalid or empty SMILES for {name!r}.")
        if len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError(f"{name!r} contains disconnected fragments; define each species separately.")
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=int(seed), useRandomCoords=True) != 0:
            raise RuntimeError(f"Could not generate a 3-D conformer for {name!r}.")
        if AllChem.MMFFHasAllMoleculeParams(mol):
            status = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        elif AllChem.UFFHasAllMoleculeParams(mol):
            status = AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
        else:
            raise ValueError(f"Neither MMFF nor UFF has all parameters for {name!r}.")
        if status != 0:
            raise RuntimeError(f"Conformer optimization did not converge for {name!r} (status {status}).")
        conf = mol.GetConformer()
        coordinates = conf.GetPositions()
        coordinates -= coordinates.mean(axis=0)
        for index, position in enumerate(coordinates):
            conf.SetAtomPosition(index, tuple(float(v) for v in position))
        self.mol_obj = mol
        self.mw_gmol = float(Descriptors.MolWt(mol))

    def copy(self):
        duplicate = object.__new__(Molecule)
        duplicate.type = self.type
        duplicate.smiles = self.smiles
        duplicate.rotate = self.rotate
        duplicate.mw_gmol = self.mw_gmol
        duplicate.mol_obj = Chem.Mol(self.mol_obj)
        return duplicate
