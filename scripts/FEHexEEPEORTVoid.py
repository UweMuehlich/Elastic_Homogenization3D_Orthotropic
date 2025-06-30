# Elastic Homogenization of a Regular Hexagonal Prism
# made of 2 different orthotropic materials
# guided by tutorial  provided by Jeremy Bleyer 
# "Periodic homogenization of linear elasticity"
# https://bleyerj.github.io/comet-fenicsx/tours/homogenization/
#                   periodic_elasticity/periodic_elasticity.html
# Ubuntu Release 18.04.6 LTS (Bionic Beaver) 64-bit
# docker version 24.0.2, build cb74dfc
#           North
#        V1 ---- V0                     
#  NWest/          \ NEast                   
#     V2            V5                       
#  SWest\          / SEast                
#        V3 ---- V4                       
#           South                      
import gmsh
import numpy as np
import dolfinx.fem.petsc
import ufl
import ast
import argparse

from dolfinx.io import XDMFFile
from dolfinx    import fem, mesh, io, default_scalar_type
from mpi4py     import MPI
from petsc4py   import PETSc

from dolfinx_mpc import MultiPointConstraint, apply_lifting
from dolfinx_mpc import assemble_matrix, assemble_vector
from dolfinx_mpc import create_sparsity_pattern

from dolfinx.mesh import meshtags

from dolfinx.io        import gmshio
from dolfinx.io.gmshio import model_to_mesh

from dolfinx_mpc import LinearProblem
#---------------------------------------------------------------------
#from   create_piecewise_constant_field import *
from   ortho_elasticity                import ortho_elasticity as oe
#---------------------------------------------------------------------
def polygon_area_2d(points):
    """
    Computes the area of a convex polygon from the coordinates of corner points
    """
    x = points[:, 0]
    y = points[:, 1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))   
#---------------------------------------------------------------------
#---------------------------------------------------------------------
if MPI.COMM_WORLD.rank == 0:
   print("===============================================")
   print(f"DOLFINx version: {dolfinx.__version__}")
   print(f"based on GIT commit: {dolfinx.git_commit_hash}")
   #print(f"of https://github.com/FEniCS/dolfinx/")
   print(f"UFL version: {ufl.__version__}")
   print("===============================================")
#---------------------------------------------------------------------
# DATA SPECIFICATION
#---------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Elastic Homogenisation HexPrism")

parser.add_argument("--ID", default="hex2matVoid", help="Problem ID")
parser.add_argument("--test" , type=int, default=-1, help="only runs case number >= 0")

args = parser.parse_args()

Case_ID, test_case = args.ID, args.test    
      
is_test   = test_case >= 0

corner_points_list = []

# read coordinates of corner points
#
with open(Case_ID+"-corners.dat", 'r') as f:
    for line in f:
        tup = ast.literal_eval(line.strip())
        corner_points_list.append(tup)

corner_points = np.array(corner_points_list)
h             = corner_points[6][2]                    # extrusion length

area = polygon_area_2d(corner_points[0:6])
print("Area of hexagon:", area)

vol = area*h 
#                                            air   fibre  lignin
# Material Constants For Materials A, B, C (inner, midddle, outer)
c_fac = 1e-12             # conversion from Pa to N / micrometer^2

#E1a,  E2a,  E3a, nu12a, nu13a, nu23a  = 10, 10,10,  0.4, 0.4, 0.4
#G12a, G13a, G23a                      = 10, 10, 10

E1b,  E2b,  E3b     = 2430080686.9651*c_fac, 2430080686.9651*c_fac, 17061238374.9446*c_fac
nu12b, nu13b, nu23b = 0.2634, 0.0364, 0.0364
G12b, G13b, G23b    = 2665000000.0*c_fac, 2665000000.0*c_fac, 1923000000.0*c_fac

E1c,  E2c,  E3c     = 1560000000.0*c_fac, 1560000000.0*c_fac, 1560000000.0*c_fac
nu12c, nu13c, nu23c = 0.3000, 0.3000, 0.3000
G12c, G13c, G23c    = 1200000000.0*c_fac, 1200000000.0*c_fac, 1200000000.0*c_fac


# Load Cases
e11, e22, e33, e12, e13, e23 = 1, 1 , 1 , 0.5, 0.5, 0.5

load_Cases = [ np.array([[e11, 0,    0], [0,   0,   0], [0,   0,0]]),
               np.array([[0  , 0,    0], [0,  e22,  0], [0,   0,0]]),
               np.array([[0  , 0,    0], [0,  0,    0], [0,   0,e33]]),
               np.array([[0  , 0,    0], [0,  0  ,e23], [0, e23,0]]),
               np.array([[0  , 0,  e13], [0,  0,    0], [e13,0,0]]),
               np.array([[0  , e12,  0], [e12,0,    0], [0,  0,0]])               
             ]

# Surface Mesh Tags for Identifying Outer Surfaces
# FrontI, FrontM, FrontO (front -inner, middle, outer-  core)
# 
NEast, North, NWest, SWest,South, SEast     = 21,16,17,18,19,20 
Front, Back                                 = 1,2
#---------------------------------------------------------------------
# READ MESH  FILE USING GMSHIO (NO PREVIOUS CONVERSION REQUIRED)
#---------------------------------------------------------------------
mesh_file = Case_ID + ".msh"      # Path to the GMSH mesh file
xdmf_file = Case_ID + ".xdmf"

domain, cell_tags, facet_tags = gmshio.read_from_msh(
              mesh_file, MPI.COMM_WORLD, 0, gdim=3)


element    = domain.geometry.cmap
ipol_order = element.degree
gdim       = domain.geometry.dim
#---------------------------------------------------------------------
# DEFINE FUNCTION SPACE
#---------------------------------------------------------------------
V = fem.functionspace(domain, ("Lagrange", ipol_order, (gdim, )))
#---------------------------------------------------------------------
# CONSTITUTIVE LAW
#---------------------------------------------------------------------
CB = oe.compliance_matrix(gdim, E1b,E2b,E3b,G12b, G13b, G23b,nu12b, nu13b, nu23b)
CC = oe.compliance_matrix(gdim, E1c,E2c,E3c,G12c, G13c, G23c,nu12c, nu13c, nu23c)

SC = oe.stiffness_matrix(CC)
SB = oe.stiffness_matrix(CB)

SB_const = oe.material_tensor(SB, domain)
SC_const = oe.material_tensor(SC, domain)

middle_vol = oe.make_cell_indicator(domain, cell_tags, tag_value = 2)
outer_vol  = oe.make_cell_indicator(domain, cell_tags, tag_value = 1)

StiffMat =  middle_vol*SB_const+outer_vol*SC_const

cell_indices = cell_tags.indices.copy() # Clone indices and tag values
tags =cell_tags.values.copy()           # 

# Create new meshtags object with modified tags
cell_tags_mod = meshtags(domain, gdim, cell_indices, tags)

#mat = create_piecewise_constant_field(
#                   domain, cell_tags_mod, {1: 1, 2: 2}, name="Material Type")

Eps  = fem.Constant(domain, np.zeros((gdim, gdim)))   # macroscopic strain
Eps_ = fem.Constant(domain, np.zeros((gdim, gdim)))   # macroscopic strain for
                                                      # post-processing
#---------------------------------------------------------------------
# VARIATIONAL PROBLEM
#---------------------------------------------------------------------
u = ufl.TrialFunction(V)                             # Trial periodic fluctuation
v = ufl.TestFunction(V)                              # Test function

# Weak form, see Jeremy Bleyer tutorial homogenization
# Function for strain and stress are defined in ortho_elasticity.py

lmbda, mu = 2*G12c*nu12c/(1-2*nu12c)*c_fac,G12c*c_fac

def epsilon(v):
    return ufl.sym(ufl.grad(v))

def sigma(v):
    eps = Eps + epsilon(v)
    return lmbda * ufl.tr(eps) * ufl.Identity(gdim) + 2 * mu * eps

LinForm = ufl.inner(oe.stress(Eps + oe.strain(u,gdim),StiffMat,gdim), oe.strain(v, gdim))
#LinForm = ufl.inner(sigma(u), epsilon(v))
a, L    = ufl.system( LinForm* ufl.dx)

#---------------------------------------------------------------------
# PERIODIC BCS with MPC from J. Dokken 
#---------------------------------------------------------------------
# fix corner points to avoid rigid body motion
#        V1 ---- V0                   
#       /          \                    
#     V2            V5                 
#       \          /                    
#        V3 ---- V4                     
#                                      
def locate_dofs(x_val):
    dofs = fem.locate_dofs_geometrical(V,lambda x: np.isclose(x[0], x_val[0]) & 
                                                   np.isclose(x[1], x_val[1]) &
                                                   np.isclose(x[2], x_val[2]))
    return dofs   

corners_coords = corner_points 
corner_dofs    = np.array([])

bcs = []

for i, coord in enumerate(corners_coords): 
    dofs = locate_dofs(coord)
    bc   = fem.dirichletbc(np.zeros((gdim,)), dofs, V)
    bcs.append(bc)

#
#         v<front>, <back>                         North
#    v1,7 ------ v0,6        y                  o ------ o       
#       /          \         |             NWest/          \NEast 
#   v2,8           v5,11     |                 o            o
#       \          /         |_______ x    SWest\          /SEast
#     v3,9 ------ v4,10                          o ------ o
#                                                  South
#  and Front - Back (vi means vertex number i)
#
# displacement vectors for mpc
s_Front_Back  = corner_points[6] - corner_points[0] 
s_South_North = corner_points[0] - corner_points[4] 
s_SWest_NEast = corner_points[5] - corner_points[3] 
s_SEast_NWest = corner_points[2] - corner_points[4] 

mpc = MultiPointConstraint(V)

periodic_relation_Front_Back =  lambda x: x + s_Front_Back[:,  np.newaxis]
periodic_relation_South_North = lambda x: x + s_South_North[:, np.newaxis]
periodic_relation_SWest_NEast = lambda x: x + s_SWest_NEast[:, np.newaxis]
periodic_relation_SEast_NWest = lambda x: x + s_SEast_NWest[:, np.newaxis]

mpc.create_periodic_constraint_topological(
    V, facet_tags, Front, periodic_relation_Front_Back, bcs)

mpc.create_periodic_constraint_topological(
    V, facet_tags, SWest, periodic_relation_SWest_NEast,bcs)

mpc.create_periodic_constraint_topological(
    V, facet_tags, SEast, periodic_relation_SEast_NWest, bcs)

mpc.create_periodic_constraint_topological(
    V, facet_tags, South, periodic_relation_South_North, bcs)

mpc.finalize()
#---------------------------------------------------------------------
# SOLVE
#---------------------------------------------------------------------
uu = fem.Function(mpc.function_space, name="Displacement")
up = fem.Function(mpc.function_space, name="Periodic_fluctuation")
ut = fem.Function(mpc.function_space, name="Total Displacements")

problem = LinearProblem(
                        a, L, mpc,
                        bcs=bcs, u=up,
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
                        )

dim_load = len(load_Cases)
CC_hom   = np.zeros((dim_load,dim_load))
#vol      = fem.assemble_scalar(fem.form(1 * ufl.dx(domain=domain)))

y = ufl.SpatialCoordinate(domain)

if is_test:
   Eps.value = load_Cases[test_case]
   up        = problem.solve()
   uu.interpolate(
       fem.Expression(
            ufl.dot(Eps, y), mpc.function_space.element.interpolation_points()
        ))
   ut.x.array[:] =  uu.x.array[:] + up.x.array[:]
   with io.VTXWriter(MPI.COMM_WORLD, "fluctuations.vtu", [up]) as vtk:
        vtk.write(0.0)  # 0.0 is the time or step

   with io.VTXWriter(MPI.COMM_WORLD, "total_disp.vtu", [ut]) as vtk:
        vtk.write(0.0)  # 0.0 is the time or step
   quit()

for nload in range(dim_load):
    Eps.value = load_Cases[nload]

    up = problem.solve()
    
    print("finished case")
    print(nload, Eps.value)
    
    for nload_ in range(dim_load):
        Eps_.value   = load_Cases[nload_]
        # Be careful here with Eps and Eps_ ! 
        stressEps    = oe.stress(Eps + oe.strain(up,gdim),StiffMat,gdim)
        #stressEps = sigma(up)
        argument     = ufl.inner(stressEps, Eps_)
        int_argument = fem.assemble_scalar(fem.form( argument * ufl.dx)) / vol

        CC_hom[nload, nload_] = int_argument
#---------------------------------------------------------------------
# RESULTS
#---------------------------------------------------------------------
#print("Aparent stiffness:")
#np.set_printoptions(precision=10, suppress=True)
#print(CC_hom)

constants = oe.extract_elastic_constants_from_stiffness(CC_hom, 1/c_fac)
print("effective constants (E,G in N/m^2")
for k, v in constants.items():
    print(f"{k}: {v:.3e}")

with io.XDMFFile(MPI.COMM_WORLD, "solution.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_meshtags(cell_tags, domain.geometry)
    xdmf.write_meshtags(facet_tags, domain.geometry)
    xdmf.write_function(up)
    xdmf.write_function(uu)
    #xdmf.write_function(mat)

with io.XDMFFile(MPI.COMM_WORLD, "mesh_with_tags.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_meshtags(facet_tags, domain.geometry)

if MPI.COMM_WORLD.rank==0:
    output_file = Case_ID + "simulation_summary.txt"

    with open(output_file, "w") as f:
        f.write("Simulation Summary\n")
        f.write("==================\n\n")

        f.write("Geometry parameters:\n")
        f.write(f" Extrusion lengh: {h}\n\n")

        f.write("Simulation parameters:\n")
        f.write(f" Interpolation order: {ipol_order} \n \n")

        
        f.write(f" Volume {vol} \n \n")

        f.write("Homogenized Elasticity Tensor (CC_hom):\n")
        for row in CC_hom:
              row_str = "  ".join(f"{val/c_fac:15.6e}" for val in row) #
              f.write(f"  {row_str}\n")
        f.write("\n")
        f.write("effective constants (E,G in N/m^2)\n")
        for k, v in constants.items():
            f.write(f"{k}: {v:.3e}\n")
        f.write("\n")
        # Get all unique facet tag values
        unique_facet_tags = np.unique(facet_tags.values)
        f.write(f"Facet tag values present: {unique_facet_tags}\n")
        f.write(f"corner points \n {corner_points} \n")   
        f.write(f"displacement vectors\n")    
        f.write(f"SFB  {s_Front_Back} \n")
        f.write(f"SSN  {s_South_North} \n")
        f.write(f"SWNE {s_SWest_NEast} \n")
        f.write(f"SENW {s_SEast_NWest} \n") 

if MPI.COMM_WORLD.rank == 0:
    print("Solution written to solution.xdmf")



