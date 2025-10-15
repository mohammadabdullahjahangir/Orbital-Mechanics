# Lambert's Problem Tools File

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import math as m
import numpy as np
import tools as t
import planetary_data as pd

pi = m.pi

# Universal Variables Method Lambert's Solver
def lamberts_universal_variables(r0, r1, deltat, tm=1, mu=pd.Earth['mu'], tol=1e-6, max_steps=200, psi=0, psi_u=4*pi**2, psi_l=-4*pi):
    # Calculate Square Root of mu Parameter
    sqrt_mu = m.sqrt(mu)
    
    # Calculate Norm of Position Vectors
    r0_norm = t.norm(r0)
    r1_norm = t.norm(r1)
    
    # Calculate Gamma Parameter
    gamma = np.dot(r0, r1) / (r0_norm * r1_norm)
    
    # Calculate Beta Parameter
    beta = tm * m.sqrt(1 - gamma ** 2)
    
    # Calculate A Parameter
    A = tm * m.sqrt(r0_norm * r1_norm * (1 + gamma))
    
    # If A = 0, Solution Can't Be Calculated
    if A == 0:
        return np.array([0, 0, 0]), np.array([0, 0, 0])  # Fixed: proper tuple syntax
    
    # Initial Values of C2 and C3 Parameters
    c2 = 0.5
    c3 = 1/6
    
    # Counter and Solved Variables
    step = 0
    solved = False
    
    # While Tolerance not Met and not at Max Steps
    for n in range(max_steps):
        # Calculate B Parameter
        B = r0_norm + r1_norm + A * (psi * c3 - 1) / m.sqrt(c2)
        
        # If A and B Parameters out of Range
        if A > 0 and B < 0:
            # Increase Lower psi Value
            psi_l += pi  # Fixed: was missing operator
            # Recalculate B
            B *= -1
        
        # Calculate Universal Variable Cubed
        chi3 = m.sqrt(B / c2) ** 3
        
        # Calculate deltat_ Variable
        deltat_ = (chi3 * c3 + A * m.sqrt(B)) / sqrt_mu  # Fixed: removed space
        
        # If Difference Between deltat Variables
        if abs(deltat - deltat_) < tol:
            # Set Solved Variable to True
            solved = True
            # Break out of Loop
            break
        
        # Check Difference between deltat and deltat_
        if deltat_ < deltat:
            # Adjust Lower psi Value
            psi_l = psi
        else:
            # Adjust Upper psi Value
            psi_u = psi
        
        # Update psi, c2 and c3 Values
        psi = (psi_u + psi_l) / 2.0  # Fixed: should be addition for bisection
        c2 = C2(psi)
        c3 = C3(psi)
        print(f"Step {n}: psi = {psi}")  # Better debug output
    
    # Check if Maximum Steps Reached
    if not solved:
        # Algorithm Didn't Converge on a psi Value
        print('Lamberts UV variables did not converge')
        return np.array([0, 0, 0]), np.array([0, 0, 0])  # Fixed: proper tuple syntax
    
    # Calculate Coefficients
    f = 1 - B / r0_norm
    g = A * m.sqrt(B / mu)
    gdot = 1 - B / r1_norm
    
    # Calculate Velocity Vectors
    v0 = (r1 - f * r0) / g  # Fixed: added space for clarity
    v1 = (gdot * r1 - r0) / g
    
    return v0, v1


# Stump Function C2
def C2(psi):
    if abs(psi) < 1e-6:  # Handle near-zero case
        return 0.5
    elif psi > 0:
        return (1 - m.cos(m.sqrt(psi))) / psi
    else:  # psi < 0
        return (m.cosh(m.sqrt(-psi)) - 1) / (-psi)


# Stump Function C3
def C3(psi):
    if abs(psi) < 1e-6:  # Handle near-zero case
        return 1/6
    elif psi > 0:
        return (m.sqrt(psi) - m.sin(m.sqrt(psi))) / (psi * m.sqrt(psi))
    else:  # psi < 0
        sqrt_neg_psi = m.sqrt(-psi)
        return (m.sinh(sqrt_neg_psi) - sqrt_neg_psi) / ((-psi) * sqrt_neg_psi)