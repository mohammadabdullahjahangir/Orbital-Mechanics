# Solar Radiation Pressure 

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
import spiceypy as spice
import spice_tools as st
 

# Total Simulation Days
tspan = 3600 * 365 * 24 # days

date0 = '2020-03-08'

h = 30e-3
w = 35e-3
A = h * w

if __name__ == '__main__':

    # Initial Conditions
    state0 = t.tle2coes('txt/geos1.txt')
    state1 = [42095.0, 0.81818, 28.5, 180.0, 298.2253, 357.857]
    
    # Null Perturbations Dictionary
    perts = null_perts()

    # Add Lunar Gravity Perturbation
    perts['SRP'] = True
    perts['A_SRP'] = A
    perts['CR'] = 1.0
    mass0 = 7000

    # Create Orbit Propagator Instance
    op0 = OP(state0, tspan, 1000, coes = True, deg = True, perts = perts, date0 = date0, mass0 = mass0, propagator = 'dopri5')
    op0.calculate_coes(parallel = False)
    op0.plot_coes(days = True, rel = True, title = 'Geostationary')
    op0.plot_3d(title = 'Geostationary Orbit')

    op1 = OP(state1, tspan, 1000, coes = True, deg = True, perts = perts, date0 = date0, mass0 = mass0, propagator = 'dopri5')
    op1.calculate_coes( parallel = False)
    op1.plot_coes(days = True, rel = True, title = 'Molniya')
    op1.plot_3d(title = 'Molniya Orbit')