# Implemented Perturbation Stop Conditions

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
tspan = 3600 * 10
dt = 100.0

cb = pd.Earth


if __name__ == '__main__':
    # Define Perturbations Dictionary
    perts = null_perts()
    perts['Aero'] = True
    perts['Cd'] = 2.2
    perts['A'] = (1e-3)**2 / 4.0
    perts['ISP'] = 4300
    perts['Thrust Direction'] = -1

    # Define Stop Conditions
    sc = {'min_alt': 300.0}

    # Initial Mass of Spacecraft
    mass0 = 50

    # Calculate Initial State Vector
    state0 = [cb['radius'] + 800, 0.0, 10.0, 0.0, 0.0, 0.0]

    op = OP(state0, tspan, dt, coes = True, deg = True, perts = perts, sc = sc, mass0 = mass0)
    op.plot_alts(show_plot = True, hours = True)
    op.plot_3d(show_plot = True)
    op.calculate_coes()
    op.plot_coes(show_plot = True, hours = True)
    op.calculate_apoapse_periapse()
    op.plot_apoapse_periapse(show_plot = True, hours = True)
