# Used tools.py file and plot_n_orbits function
import sys 
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
import tools as t
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP


tspan = 24 * 60 * 60  
dt = 100.0

if __name__ == '__main__':
    r_mag = pd.Earth['radius'] + 500.0
    v_mag = np.sqrt(pd.Earth['mu'] / r_mag)

    r0 = np.array([r_mag, 0.0, 0.0])
    v0 = np.array([0.0, v_mag, 0.0])

    r_mag = pd.Earth['radius'] + 1000.0
    v_mag = np.sqrt(pd.Earth['mu'] / r_mag) * 1.3 

    r00 = np.array([r_mag, 0.0, 0.0])
    v00 = np.array([0.0, v_mag, 0.3])

    op0 = OP(r0, v0, tspan, dt)
    op00 = OP(r00, v00, tspan, dt)
    op0.propagate_orbit()
    op00.propagate_orbit()  
    t.plot_n_orbits([op0.rs, op00.rs], ['Orbit 1', 'Orbit 2'])