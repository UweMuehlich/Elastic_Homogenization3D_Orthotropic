import numpy as np

def extract_orthotropic_constants(C):
    """
    Extract orthotropic engineering constants from the Stiffness matrix C 
    in Voigt notation.
    
    Parameters:
        C: numpy.ndarray of shape (6, 6), the Stiffness matrix in Voigt notation.
        
    Returns:
        A dictionary with E1, E2, E3, nu12, nu13, nu23, G12, G13, G23
    """
    if C.shape != (6, 6):
        raise ValueError("Stiffness matrix must be 6x6 in Voigt notation.")

    # Compute compliance matrix
    S = np.linalg.inv(C)

    # Engineering constants
    E1 = 1.0 / S[0, 0]
    E2 = 1.0 / S[1, 1]
    E3 = 1.0 / S[2, 2]

    nu12 = -S[0, 1] * E1
    nu13 = -S[0, 2] * E1
    nu23 = -S[1, 2] * E2

    G12 = 1.0 / S[3, 3]
    G13 = 1.0 / S[4, 4]
    G23 = 1.0 / S[5, 5]

    return {
        'E1': E1,
        'E2': E2,
        'E3': E3,
        'nu12': nu12,
        'nu13': nu13,
        'nu23': nu23,
        'G12': G12,
        'G13': G13,
        'G23': G23
    }

# Example usage:
if __name__ == "__main__":
    # Example orthotropic stiffness matrix (units: GPa)
    S_example = np.array([
        [150,  5,  5,   0,   0,   0],
        [  5, 80,  5,   0,   0,   0],
        [  5,  5, 50,   0,   0,   0],
        [  0,  0,  0, 30,   0,   0],
        [  0,  0,  0,   0, 20,   0],
        [  0,  0,  0,   0,   0, 10]
    ])
    
    S_fibre =  np.array([
    [2.654e9,  7.305e8, 8.656e8,     0.0,     0.0,     0.0],
    [7.305e8,  2.654e9, 8.656e8,     0.0,     0.0,     0.0],
    [8.656e8,  8.656e8, 1.7504e10,   0.0,     0.0,     0.0],
    [0.0,      0.0,     0.0,         2.665e9, 0.0,     0.0],
    [0.0,      0.0,     0.0,         0.0,     2.665e9, 0.0],
    [0.0,      0.0,     0.0,         0.0,     0.0,     1.923e9]
     ])

    S_lignin = np.array([
      [2.1e9, 9e8,   9e8,   0.0,   0.0,   0.0],
      [9e8,   2.1e9, 9e8,   0.0,   0.0,   0.0],
      [9e8,   9e8,   2.1e9, 0.0,   0.0,   0.0],
       [0.0,   0.0,   0.0,   1.2e9, 0.0,   0.0],
       [0.0,   0.0,   0.0,   0.0,   1.2e9, 0.0],
       [0.0,   0.0,   0.0,   0.0,   0.0,   1.2e9]
    ])
    

    constantsf = extract_orthotropic_constants(S_fibre)
    print("FIBRE")
    for name, value in constantsf.items():
        print(f"{name}: {value:.4f}")

    constantsl = extract_orthotropic_constants(S_lignin)
    print("LIGNIN")
    for name, value in constantsl.items():
        print(f"{name}: {value:.4f}")
