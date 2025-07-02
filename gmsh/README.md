# Hexagonal Prism Mesh Generator for FEniCSx

This part of the repository provides a **parametric Gmsh mesh generator** for a hexagonal prism unit cell with up to three concentric regions (outer ring, middle ring, inner core). It uses **Gmsh’s Python API** and converts the mesh to XDMF for use with **FEniCSx**.

## 📌 Main Script

- **Script:** `hexWoodVoidGmshAPI.py`
- **Purpose:**  
  - Generates a 3D hexagonal prism geometry with specified dimensions.
  - Supports creating outer, middle, and optional inner rings.
  - Tags physical groups for volumes, faces, and edges for applying boundary conditions in FEniCSx.
  - Exports mesh in `.msh` and `.xdmf` formats.
  - Saves corner coordinates for reference.
  - Helper functions which do not depend on gmsh are collected in ./helpers/parametric_hexagon.py

---

<img src="../docs/images/hex_cell.svg" alt="Hexagonal Prism Unit Cell" width="400"/>

- **Top view:** three concentric hexagonal regions: Outer Ring, Middle Ring, and Inner Core (optional).
- **Side view:** the extrusion height specified by `--height`.

---


### ⚙️ Command-line Options 

| Option         | Description                                        | Default               |
|----------------|----------------------------------------------------|-----------------------|
| `--ID`         | Problem ID used for filenames                      | `three_hex_rings`     |
| `--T`          | Width of hexagon                                   | `2.0`                 |
| `--D`          | Length of slanted edge                             | `2.5`                 |
| `--angle`      | Hexagon angle in degrees                           | `60`                  |
| `--thick`      | Thickness of the middle ring                       | `0.4`                 |
| `--thickOut`   | Thickness of the outer ring                        | `0.2`                 |
| `--height`     | Extrusion length of the prism                      | `1.0`                 |
| `--MeshSize`   | Global mesh size                                   | `0.5`                 |
| `--IntOrder`   | Finite element order for mesh generation (1 or 2)  | `1`                   |
| `--inner`      | Include the inner core                             | `False` (flag option) |
| `--prisms`     | Use prism elements instead of tetrahedra           | `False` (flag option) |


### 📂 Output Files

All output files are saved in the `mesh/` folder:

| File                      | Description                                                      |
|---------------------------|------------------------------------------------------------------|
| `<ID>.msh`                | Gmsh mesh file with physical groups (for inspection in Gmsh GUI) |
| `<ID>.xdmf`               | Converted XDMF mesh for use with FEniCSx simulations             |
| `<ID>-corners.dat`        | Text file with coordinates of the hexagon corner points (front and back faces) |


### 1️⃣ Install dependencies

```bash
pip install gmsh meshio

