# Changelog

## 2.0.0

- Incremental Sobol insertion with largest-first placement and periodic vdW clash rejection.
- Targeted density and increasing-density search, saving the best accepted configuration.
- Independent multicore placement attempts, stopping when a successful result arrives.
- Seeded packing, actual-density reports, and explicit failure reporting.
- Installable Python package, editable input template, and command-line interface.
- Legacy input-name translation; acceptance now uses geometric distances instead of pseudo-LJ energy.

The purpose remains SMILES-to-PDB initial system generation. Old Python internals
are not a supported API, and the packing criterion and random-seed trajectories
have changed. Existing input files can be loaded through `--config`; see the README.
