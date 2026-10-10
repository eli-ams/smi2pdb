"""User configuration for smi2pdb 2.0.0; lengths are in angstroms."""
output_dir = "output/"  # Relative to the directory containing your input file.
mixtures_needed = 1
run_mode = "targeted"  # "targeted" or "density_search".
target_density_kg_m3 = 1000.0
initial_density_fraction = 0.30
density_search_increment = 0.10  # Search starts at initial_density_fraction.
density_search_max_fraction = 1.0  # Search never exceeds 1.0 * target density.
box_shape = "cubic"  # Or "elongated_z" for scymol's box shape.
random_seed = 42
number_of_trials = 10  # Full packing attempts per box size, per mixture.
packing_workers = 0  # 0: automatic (CPUs minus one); 1: serial; or e.g. 4.
candidates_per_molecule = 32  # Shared budget: max(512, molecule_count * this).
clash_scale = 0.72  # Minimum intermolecular distance / sum of vdW radii.
boundary_padding = 0.0  # Extra margin beyond molecular vdW extent.
max_box_expansions = 0  # Strict density by default; set to 5 for fallback.
box_expansion_factor = 1.20  # Expands all three dimensions; lowers density.
molecules_dictionary = {
    "phenolic_asphaltene": {
        "nbr_of_mols": 500,
        "smiles": "CC(C)Cc1cc2c3c(cc4c(CC(CC4CC)CCCC)c3c1)c1cc(O)cc3c1c2cc(CCC)c3CCC(C)C",
        "rotate": True,
    },
}
AVOGADRO_NUMBER = 6.02214076e23
