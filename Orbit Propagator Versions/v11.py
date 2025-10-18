# Modelling Spacecrafts with Aerodynamic Drag

import sys 
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

# Time Parameters
tspan = 24 * 3600
dt = 100.0

cb = pd.Earth

if __name__ == '__main__':
    # Define Perturbations Dictionary
    perts = null_perts()
    perts['Aero'] = True
    perts['Cd'] = 2.2
    perts['A'] = (1e-3)**2 / 4.0

    # Initial Mass of Spacecraft
    mass0 = 10

    # Perigee and Apogee Altitudes
    rp = 215 + cb['radius']
    ra = 300 + cb['radius']

    # Orbital Elements
    raan = 340.0
    i = 65.2
    aop = 58.0
    ta = 332.0

    # Calculate Other Orbital Elements
    a = (ra + rp) / 2.0
    e = (ra - rp) / (ra + rp)

    # Initial State Vector
    state0 = [a, e, i, ta, aop, raan]

    op = OP(state0, None, None, tspan, dt, coes = True, deg = True, perts = perts)
    op.plot_alts(show_plot = True, hours = True)
    op.plot_3d(show_plot = True)
    op.calculate_coes(degrees = True)
    op.plot_coes(show_plot = True, hours = True)
    op.calculate_apoapse_periapse()
    op.plot_apoapse_periapse(show_plot = True, hours = True)