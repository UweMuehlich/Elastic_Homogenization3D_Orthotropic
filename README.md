# Effective Elastic Properties of Wood at Micrometer Scale

This repository contains a FEniCSx-based simulation for determining the effective elastic properties of wood at the micrometer scale. To this end, a unit cell geometry in the form of a hexagonal prism is considered. The cell features a prescribed wall thickness, with an outer  and an inner layer.

The main challenges of this simulation are:

- Modeling the three-dimensional geometry with periodic boundary conditions
- Combination of different materials, including orthotropic materials 
- Capturing the effective  orthotropic material behavior

<img src="./docs/images/total_disp_case_5.png"
     alt="Hexagonal Prism Unit Cell"
     width="250"/>


## Project structure

- `data/`   : Example data for material properties
- `scripts/`: Python simulation scripts for FEniCSx 
- `gmsh/`   : Python scripts for mesh generation using gmsh API
- `results/`: Output files (add to .gitignore)
- `docker/` : Docker instructions
- `docs/`   : Documentation files
- `tests/`  : Unit or integration tests


```
