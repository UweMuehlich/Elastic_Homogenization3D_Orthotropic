# Elastic Homogenization of a Regular Hexagonal Prism

The Python script FEHexEEPEORTVoid.py runs a FEniCSx-based simulation for determining the effective elastic properties of a hexagonal unit cell composed of two orthotropic materials. The structure is representative of biological or composite microstructures with layered or porous architecture. Periodic boundary conditions are imposed using multi-point constraints (MPC).

## Overview

The simulation computes the homogenized (effective) elasticity tensor for a hexagonal prism by solving a sequence of linear elasticity problems under six independent strain states. The method is based on periodic homogenization theory and builds upon the tutorial by [Jeremy Bleyer](https://bleyerj.github.io/comet-fenicsx/tours/homogenization/periodic_elasticity/periodic_elasticity.html).

---

## Requirements

- Python ≥ 3.10
- FEniCSx (compatible version with `dolfinx_mpc`)
- GMSH ≥ 4.11 (for mesh generation)
- MPI (e.g. `mpich` or `openmpi`)
- `petsc4py`, `numpy`, `ufl`

---

## Input Files

- `hex2matVoid.msh`: GMSH-generated mesh file with facet and volume tags.
- `hex2matVoid-corners.dat`: List of 3D corner point coordinates for the unit cell.

---

## Command-line Arguments

The script uses `argparse` to define two command-line arguments:

| Argument  | Description                                                                                          | Default         |
|-----------|------------------------------------------------------------------------------------------------------|-----------------|
| `--ID`    | Problem ID used to locate the input files (`.msh` mesh and `-corners.dat` corners). Example: `hex2matVoid` means the mesh should be `hex2matVoid.msh`. | `"hex2matVoid"` |
| `--test`  | Run only a single test case (`0`–`5` for the 6 independent load cases). If `--test` ≥ 0, only that case is run and the fluctuation and total displacement fields are saved to `.vtu` files. | `-1` (run all cases) |

---
## Output Files

After running the simulation, the following files will be generated:

- `solution.xdmf`: Complete mesh and displacement fields (for Paraview visualization).
- `fluctuations.vtu`: Periodic fluctuation displacement field (only if `--test` ≥ 0).
- `total_disp.vtu`: Total displacement field = macro + fluctuation (only if `--test` ≥ 0).
- `hex2matVoidsimulation_summary.txt`:  
  Summary file containing:
  - Geometry parameters (extrusion height, unit cell volume)
  - Interpolation order
  - Homogenized elasticity tensor (`C^hom`)
  - Extracted elastic constants (converted back to physical units)
  - Facet tag values used for periodic constraints
  - Coordinates of corner points and periodic displacement vectors

The summary file is named using the `--ID` argument, followed by `simulation_summary.txt`.  
For example, if `--ID hex2matVoid`, the file will be:
```text
hex2matVoidsimulation_summary.txt

 ---
## Example Usage

Run the full homogenization (all load cases):

```bash
python3 FEHexEEPEORTVoid.py --ID hex2matVoid

Run only test case 0 (e11 loading):

```bash
python3 FEHexEEPEORTVoid.py --ID hex2matVoid --test 0

## Helper Functions and Modules

The main script relies on a custom module `ortho_elasticity.py`, which provides helper functions for:

- Defining the orthotropic compliance and stiffness matrices
- Computing stress and strain tensors for the variational form
- Extracting effective elastic constants from the homogenized stiffness tensor
- Creating piecewise constant material tensors and cell indicators

These helper functions are located in the `ortho_elasticity` directory or file.

➡️ **Note:** A separate `README.md` is provided in the `ortho_elasticity/` folder. Please refer to it for details on the available functions, their usage, and examples.

---




