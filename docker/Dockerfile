# FEniCSx Dockerfile
# Run with: docker build -t dolfinx_mpc_gmsh -f docker/Dockerfile .

FROM ghcr.io/jorgensd/dolfinx_mpc:v0.9.1

# Install GMSH with Python API

RUN apt-get update && apt-get install -y gmsh \
    && pip3 install gmsh

# Add any extra Python libraries you need
# RUN pip3 install pyvista matplotlib

# Set working directory inside container
WORKDIR /root/fenics

# (Optional) default command to run bash
CMD ["bash"]

