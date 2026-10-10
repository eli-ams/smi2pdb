# smi2pdb

**From SMILES strings to packed molecular systems.**

smi2pdb is a Python tool for preparing molecular configurations for atomistic
simulations. Define your molecules, their counts, and the desired density;
smi2pdb generates three-dimensional conformers, packs them into a periodic box,
and exports the complete system as a PDB file.

Version **2.0.0** introduces incremental steric clash rejection, density search,
and multicore packing for efficient preparation of organic mixtures.

![smi2pdb 2.0.0 workflow: molecular definitions, conformer generation, density selection, incremental Sobol packing, parallel attempts, optional density search, and PDB export](https://raw.githubusercontent.com/eli-ams/smi2pdb/master/readme_resources/workflow_2_0.png)

*The smi2pdb 2.0.0 workflow. Molecular illustrations and density-search values are schematic.*

## Highlights

- **SMILES to 3D:** generate explicit-hydrogen conformers with RDKit and optimize them with MMFF or UFF.
- **Incremental packing:** place larger molecules first, retaining successful insertions while retrying candidates that clash.
- **Two density modes:** pack at a chosen density or search progressively higher densities up to a specified fraction of the reference density.
- **Multicore attempts:** run independent placements in parallel and retain the first complete packing.
- **Traceable output:** save coordinates, input settings, seeds, achieved densities, and packing statistics.

## Installation

Requires **Python 3.10 or newer**. NumPy and RDKit are installed as dependencies.

```sh
python -m pip install smi2pdb
```

To install from a downloaded or cloned source checkout, run this in its root directory:

```sh
python -m pip install .
```

## Quick start

Create an editable configuration:

```sh
smi2pdb --init inputs.py
```

Edit the molecular definitions and packing settings in `inputs.py`, then run:

```sh
smi2pdb --config inputs.py
```

The source checkout includes an `inputs.py` example; use it directly or create a
separate configuration. `--init` preserves existing files. The same commands are
available through `python -m smi2pdb`.

### Define a mixture

For example, save the following as `inputs.py`:

```python
output_dir = "output/"
mixtures_needed = 1

run_mode = "targeted"
target_density_kg_m3 = 1000.0
initial_density_fraction = 0.30
packing_workers = 0

molecules_dictionary = {
    "hexane": {
        "smiles": "CCCCCC",
        "nbr_of_mols": 200,
        "rotate": True,
    },
    "toluene": {
        "smiles": "Cc1ccccc1",
        "nbr_of_mols": 100,
        "rotate": True,
    },
}
```

This requests 200 hexane molecules and 100 toluene molecules at **300 kg/m³**:
30% of the specified 1000 kg/m³ reference density. Settings omitted from the file
use the package defaults. Each species needs a connected SMILES definition,
a positive molecule count, and a rotation setting.

Input files are Python configuration files. Relative output paths resolve beside
the configuration file, so each project can keep its inputs and results together.

## Choose a run mode

### Targeted density

Use `targeted` to pack at a specific fraction of your reference density:

```python
run_mode = "targeted"
target_density_kg_m3 = 1000.0
initial_density_fraction = 0.50
```

Here, the requested packing density is **500 kg/m³**. The box stays at the requested
size by default. Optional box expansion can be enabled with `max_box_expansions`;
the resulting density is always recorded in the packing report.

### Increasing-density search

Use `density_search` to find the highest density successfully packed within your
chosen range and attempt budget:

```python
run_mode = "density_search"
target_density_kg_m3 = 1000.0
initial_density_fraction = 0.30
density_search_increment = 0.10
density_search_max_fraction = 1.0
```

The search tries **0.30ρ → 0.40ρ → 0.50ρ → … → 1.00ρ**, where ρ is
`target_density_kg_m3`. It stops when a density exhausts its packing attempts or
when the upper endpoint is reached. Each step starts from the same molecular
conformers with new placement seeds.

Every successful structure is saved. If 0.50ρ succeeds and 0.60ρ fails, the main
PDB retains the 0.50ρ configuration. Search mode disables box expansion and keeps
the same clash criterion throughout. `mixtures_needed` selects how many independent
searches to perform.

## Multicore packing

```python
packing_workers = 0  # Automatic
```

Automatic mode uses the available CPU count minus one, with at least one worker,
and is capped by the number of permitted attempts. Set `packing_workers = 4` for
up to four workers or `packing_workers = 1` for serial operation.

Workers try independent arrangements at the same density. The first complete
packing wins, and remaining workers stop. `number_of_trials` is the **total**
attempt budget per box size, not a per-worker budget. Density steps remain sequential.

You can override mode and worker count without editing the input file:

```sh
smi2pdb --config inputs.py --mode density_search --workers 4
```

Use one worker when deterministic winner selection is important. Parallel selection
can depend on scheduling; reports record the winning attempt seed. Each worker
holds its own packing data, so worker count also controls memory use.

## How packing works

Molecules are inserted from largest spatial extent to smallest. For each insertion,
smi2pdb proposes a Sobol position and, when enabled, a uniformly sampled rotation.
It accepts the placement only when all candidate atoms satisfy the configured
separation from previously placed molecules:

**dᵢⱼ ≥ clash_scale × (rᵢ + rⱼ)**

Here, rᵢ and rⱼ are the atoms' van der Waals radii. The default `clash_scale` is
**0.72**. Periodic minimum-image distances and spatial bins make these local
checks efficient. Rejected candidates are replaced while accepted placements
remain in the box.

Packing produces initial coordinates for subsequent force-field assignment,
minimization, and equilibration. The density-search result is the highest density
found under the selected geometry, clash criterion, and search budget.

![Example molecular system visualized in OVITO](https://raw.githubusercontent.com/eli-ams/smi2pdb/master/readme_resources/mixture_ortho.png)

*Example molecular system from the project's original illustrations, visualized in OVITO.*

## Main settings

| Setting | Purpose | Default |
| --- | --- | --- |
| `run_mode` | `targeted` or `density_search` | `"targeted"` |
| `target_density_kg_m3` | Reference density, ρ, in kg/m³ | `1000.0` |
| `initial_density_fraction` | Targeted fraction or search starting point | `0.30` |
| `density_search_increment` | Fraction added at each search step | `0.10` |
| `density_search_max_fraction` | Search endpoint, at most 1.0 | `1.0` |
| `mixtures_needed` | Number of mixtures or independent searches | `1` |
| `packing_workers` | Automatic, serial, or explicit worker count | `0` |
| `number_of_trials` | Full packing attempts per box size | `10` |
| `candidates_per_molecule` | Shared candidate budget multiplier | `32` |
| `clash_scale` | Minimum-distance factor applied to vdW radius sums | `0.72` |
| `box_shape` | `cubic` or `elongated_z` | `"cubic"` |
| `max_box_expansions` | Optional expansion rounds in targeted mode | `0` |
| `random_seed` | Base seed for conformers and packing attempts | `42` |

## Output

Each job creates a new numbered directory:

```text
output/1/
├── settings.json
├── log.txt
├── job_summary.json
└── mixture_1/
    ├── packing.json
    └── pdbs/
        └── mixture_1.pdb
```

Density searches also write `density_search.json` and a `density_steps` directory
containing the individual results. The main PDB and packing report correspond
to the highest successful step.

PDB files contain cell dimensions, atomic coordinates, and bond connectivity.
Open them in molecular visualization tools such as OVITO, VMD, or PyMOL.
Reports include requested and actual density, worker settings, attempt seeds,
and rejection statistics. Only complete packings are exported. Standard PDB
limits of 99,999 atoms and 9,999 molecular residues apply.

## Migrating from 1.0

Load an existing input file with `smi2pdb --config path/to/inputs.py`.
Legacy `density` or `initial_density` values are interpreted as absolute densities
in kg/m³ when `target_density_kg_m3` is absent. Older count and trial names are
translated automatically.

The packing criterion now uses inter-molecular vdW distances. Legacy
`energy_threshold` and `cutoff_distance` values are ignored with a warning;
`layer_offset` is translated to additional boundary padding. The new molecular
extent margins also apply.

## Development

```sh
python -m pip install -e .
python -m unittest discover -s tests
```

With dependencies available, `python run_smi2pdb.py --config inputs.py` runs directly
from a source checkout.

## Citation and credits

If you use smi2pdb in research, please cite:

> SMI2PDB: A self-contained Python tool to generate atomistic systems of organic
> molecules using their SMILES notations.

Created by **Eli I. Assaf, Xueyan Liu, Peng Lin, and Sandra Erkens**, TU Delft.
Version 2.0.0 incorporates packing ideas developed in the author's scymol project.

Licensed under the [GNU General Public License v3.0](https://github.com/eli-ams/smi2pdb/blob/master/LICENSE).
