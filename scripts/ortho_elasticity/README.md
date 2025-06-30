# Orthotropic Elasticity Utilities for FEniCSx

This repository contains a Python module providing reusable functions for modeling **orthotropic elasticity** problems in **2D** and **3D** using [**FEniCSx**](https://fenicsproject.org/).

---

## 📄 Overview

`ortho_elasticity.py` includes:

- ✅ Computation of **compliance** and **stiffness** matrices for orthotropic materials
- ✅ Conversion of stiffness matrices to `dolfinx.fem.Constant` objects
- ✅ Functions for computing **strain** and **stress** in tensorial and Voigt notation
- ✅ Support for **subdomain-specific** material properties via **indicator functions**
- ✅ Utilities to **extract engineering constants** from a given stiffness matrix

This makes it suitable for **finite element simulations** of materials like **wood**, **composites**, or **layered structures**.

---

## 🧩 Main Functions

| Function | Description |
| -------- | ----------- |
| `compliance_matrix(dim, E1, E2, ...)` | Returns the compliance matrix (Voigt notation) for orthotropic materials in 2D or 3D |
| `stiffness_matrix(S)` | Computes the stiffness matrix by inverting the compliance matrix |
| `material_tensor(C, mesh)` | Converts a NumPy stiffness matrix to a `dolfinx.fem.Constant` |
| `strain(u, dim)` | Computes the symmetric strain tensor from the displacement field |
| `stress(eps, C, dim)` | Computes the stress tensor from strain using the constitutive law |
| `make_cell_indicator(mesh, cell_tags, tag_value)` | Creates an indicator function for tagged subdomains |
| `extract_elastic_constants_from_stiffness(C)` | Extracts engineering constants (E, ν, G) from a stiffness matrix |

---


## ⚙️ Requirements

- Python ≥ 3.8
- [NumPy](https://numpy.org/)
- [UFL](https://fenics.readthedocs.io/projects/ufl/en/latest/)
- [FEniCSx](https://fenicsproject.org/)
- [petsc4py](https://petsc4py.readthedocs.io/en/stable/)

Install FEniCSx according to its [official installation guide](https://docs.fenicsproject.org/dolfinx/latest/installation.html).

---

## License

This project is licensed under the [MIT License](LICENSE).

Authored by **Uwe Mühlich**, 2025.

## Acknowledgments
This module was developed  with drafting help from OpenAI’s ChatGPT.


