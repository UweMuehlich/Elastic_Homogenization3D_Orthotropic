"""
Orthotropic Elasticity Utilities for Dolfinx
============================================

This module provides reusable functions for solving orthotropic elasticity problems 
in 2D or 3D using Dolfinx. It includes:

Functions:
----------
- compliance_matrix: Compute the compliance matrix from engineering constants.
- stiffness_matrix:  Compute the stiffness matrix by inverting the compliance matrix.
- material_tensor:   Convert stiffness matrix to a fem.Constant for variational forms.
- strain: Compute symmetric strain tensor (tensorial or vectorial form).
- stress: Compute stress vector or tensor from strain using stiffness tensor.

Features:
---------
- Supports both 2D (plane stress/strain) and full 3D orthotropic elasticity.
- Uses Voigt notation for strain and stress vectors.
- Compatible with subdomain-specific constitutive laws using indicator functions.

Intended for:
-------------
- Finite element simulation of orthotropic materials (e.g., wood, composites).
- Educational or research use in elasticity modeling with FEniCSx.

Author: Uwe Muhlich, ChatGPT 
Date: 31.05.2025
"""

import numpy as np
import ufl

from dolfinx  import fem
from petsc4py import PETSc

def is_vectorial(expr):
    """Return True if expr is a vector representing strain or stress in 
       Voigt notation."""
    return expr.ufl_shape in [(3,), (6,)]

def is_tensorial(expr):
    """Return True if expr is a second-order tensor (2×2 or 3×3 matrix)."""
    return expr.ufl_shape in [(2, 2), (3, 3)]

def strain_as_vector(eps):
    """Return strain tensor object as vector in Voigt notation""" 
    if is_tensorial(eps):
       if   eps.ufl_shape in [(2,2)]:
            return ufl.as_vector([eps[0, 0], eps[1, 1], 2 * eps[0, 1]])
       elif eps.ufl_shape in [(3,3)]:
            return ufl.as_vector([eps[0, 0], eps[1, 1], eps[2, 2], 
                                  2*eps[1, 2], 2*eps[0, 2], 2*eps[0, 1]])
    return(eps)

def compliance_matrix(dim, E1, E2, E3=None, G12=None, G13=None, G23=None, 
                      nu12=None, nu13=None, nu23=None):
    """
    Returns the compliance matrix (Voigt notation) for orthotropic material.
    """
    if dim == 2:
        S = np.array([
            [1. / E1, -nu12 / E1,        0.],
            [-nu12 / E1, 1. / E2,        0.],
            [0.,        0.,        1. / G12]
        ])
    elif dim == 3:
        # symmetry: nu21 = (E2/E1)*nu12, etc., assumed for inversion consistency
        nu21 = (E2 / E1) * nu12
        nu31 = (E3 / E1) * nu13
        nu32 = (E3 / E2) * nu23
        S = np.array([
            [1./E1, -nu21/E2, -nu31/E3,      0.,     0.,     0.],
            [-nu12/E1, 1./E2, -nu32/E3,      0.,     0.,     0.],
            [-nu13/E1, -nu23/E2, 1./E3,      0.,     0.,     0.],
            [0.,      0.,      0.,      1./G23,     0.,     0.],
            [0.,      0.,      0.,      0.,     1./G13,     0.],
            [0.,      0.,      0.,      0.,     0.,     1./G12],
        ])
    else:
        raise ValueError("Only dimensions 2 or 3 are supported.")
    return S

def stiffness_matrix(S):
    """
    Inverts the compliance matrix to compute the stiffness matrix.
    """
    C = np.linalg.inv(S)
    return C

def material_tensor(C, mesh):
    """
    Converts numpy stiffness matrix to a fem.Constant usable in variational forms.
    """
    flat = tuple(tuple(C[i, j] for j in range(C.shape[1])) for i in range(C.shape[0]))
    return fem.Constant(mesh, PETSc.ScalarType(flat))

def extract_elastic_constants_from_stiffness(C, conversion_factor=1.0):
    """
    Extract elastic constants (E1, E2, E3, nu12, nu13, nu23, G12, G13, G23)
    from a given 6x6 stiffness matrix in Voigt notation.

    Parameters
    ----------
    C : (6,6) ndarray
        Stiffness matrix in Voigt notation (units: stress/strain)

    conversion_factor : float, optional
        Factor to multiply Young's moduli (E*) and shear moduli (G*)
        to convert units (e.g., from Pa to GPa or to simulation units)

    Returns
    -------
    elastic_constants : dict
        Dictionary containing:
        E1, E2, E3         - Young's moduli
        nu12, nu13, nu23   - Poisson's ratios (engineering)
        G12, G13, G23      - Shear moduli
    """
    S = np.linalg.inv(C)

    # Young's moduli
    E1 = conversion_factor / S[0, 0]
    E2 = conversion_factor / S[1, 1]
    E3 = conversion_factor / S[2, 2]

    # Poisson's ratios
    nu12 = -S[0, 1] / S[0, 0]
    nu13 = -S[0, 2] / S[0, 0]
    nu21 = -S[1, 0] / S[1, 1]
    nu23 = -S[1, 2] / S[1, 1]
    nu31 = -S[2, 0] / S[2, 2]
    nu32 = -S[2, 1] / S[2, 2]

    # Shear moduli
    G12 = conversion_factor / S[3, 3]
    G13 = conversion_factor / S[4, 4]
    G23 = conversion_factor / S[5, 5]

    return {
        'E1': E1, 'E2': E2, 'E3': E3,
        'nu12': nu12, 'nu13': nu13, 'nu23': nu23,
        'nu21': nu21, 'nu31': nu31, 'nu32': nu32,
        'G12': G12, 'G13': G13, 'G23': G23
    }


def strain(u, dim):
    eps = ufl.sym(ufl.grad(u))
    return eps

def stress(eps, C, dim=2):
    """
    Returns the stress tensor using the constitutive law 
    stress = Compliance : strain in Voigt notation
    """
    if is_tensorial(eps): epsilon = strain_as_vector(eps)
    sigma_vec = ufl.dot(C, epsilon)

    if   dim == 2:
            return ufl.as_matrix([[sigma_vec[0], sigma_vec[2]],
                                  [sigma_vec[2], sigma_vec[1]]])
    elif dim == 3:
            return ufl.as_matrix([
                [sigma_vec[0], sigma_vec[5], sigma_vec[4]],
                [sigma_vec[5], sigma_vec[1], sigma_vec[3]],
                [sigma_vec[4], sigma_vec[3], sigma_vec[2]]
            ])


def make_cell_indicator(my_mesh, cell_tags, tag_value):
    """
    Create a DG0 scalar indicator function that is 1 in cells with tag_value,
                                                   0 elsewhere.
    """
    # Create DG0 element and function space
    V0 = fem.functionspace(my_mesh, ("DG", 0))
    
    # Allocate array and set indicator = 1 where tag matches
    values = np.zeros(my_mesh.topology.index_map(my_mesh.topology.dim).size_local,
                                                 dtype=np.float64)
    values[cell_tags.find(tag_value)] = 1.0
    
    # Assign to fem.Function
    indicator = fem.Function(V0)
    indicator.x.array[:] = values
  
    return indicator

