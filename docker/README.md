# 🐳 Reproducible Docker environment for FEniCSx + MPC + GMSH

This Docker image provides:
- FEniCSx with multi-point constraints (`dolfinx_mpc`)
- GMSH with Python API (pinned version)
- Other pinned Python packages (e.g., PyVista)

This ensures that your mesh generation and finite element simulations run in
a **fully reproducible** environment, regardless of your local OS.

---

## ✅ How to build the Docker image (optional)

You can build the image locally if you want to test changes:

```bash
docker build -t dolfinx_mpc_gmsh -f docker/Dockerfile .

# Make sure the script fenicsx_mpc_gmsh.sh is executable:
chmod +x docker/fenicsx_mpc_gmsh.sh

# Run the container with your project mounted:
./docker/fenicsx_mpc_gmsh.sh


