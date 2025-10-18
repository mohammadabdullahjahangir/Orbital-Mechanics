# n-Body Perturbations

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


# Central Body
cb = pd.Earth

# Total Simulation Time and Steps
tspan = 3600.0 * 24.0 * 100.0
dt = 5000

date0 = '2020-02-23'

if __name__ == '__main__':
    # Initial Conditions

    iss_coes = t.tle2coes('txt/iss.txt')

    state0 = [42164.0, 0.001, 0.0, 0.0, 0.0, 0.0]

    # Null Perturbation Theory
    perts = null_perts()

    #Add Lunar Gravity Perturbation
    perts['N_Bodies'] = [pd.Moon]

    # Create Orbit Propagator Object for Geostationary and ISS Orbit
    op0 = OP(state0, tspan, dt, coes = True, date0 = date0, perts = perts, propagator = 'dopri5')
    op_ISS = OP(iss_coes, tspan, 1000, coes = True, date0 = date0, deg = True, perts = perts, propagator = 'dopri5')

    op_ISS.calculate_coes(parallel = False)
    op_ISS.plot_coes(days = True, show_plot = True)

    op0.calculate_coes(parallel = False)
    op0.plot_coes(days = True, show_plot = True)

    t.plot_n_orbits([op0.rs, op_ISS.rs, op0.perts['N_Bodies'][0]['states'][:,:3]], ['Geostationary', 'ISS', 'Moon'], show_plot = True, ER = True, cb = pd.Earth, figsize = (10,10))