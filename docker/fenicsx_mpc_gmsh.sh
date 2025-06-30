#!/bin/bash

exec docker run -v $(pwd):/root/fenics/shared -ti dolfinx_mpc_gmsh
