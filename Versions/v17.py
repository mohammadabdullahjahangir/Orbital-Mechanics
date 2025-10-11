# Orbit Eccentricity and Escape Velocities

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

# Time Parameters
tspan = 3600 * 24 * 3 
dt = 100

date0 = '2020-04-03'

if __name__ == '__main__':

    # Calculate Initial State Vectors
    r0_apogee, v0_apogee = t.coes2rv([cb['radius'] + 26600, 0.74, 35.0, 180.0, 0.0, 0.0])
    r0_perigee, v0_perigee = t.coes2rv([cb['radius'] + 26600, 0.74, 35.0, 180.0, 0.0, 0.0])
    coes0_rotated = [cb['radius'] +26600, 0.74, 35.0, 0.0, 45.0, 0.0]

    # Calculate Circular Velocities at Initial Positions
    v_circ_apogee = (pd.Earth['mu'] / t.norm(r0_apogee)) ** 0.5
    v_circ_perigee = (pd.Earth['mu'] / t.norm(r0_perigee)) ** 0.5

    # Calculate Normed Velocity Vectors
    v0_apogee_normed = t.normed(v0_apogee)
    v0_perigee_normed = t.normed(v0_perigee)

    # Calculate Escape Velocities
    esc_v_apogee = t.esc_v(t.norm(r0_apogee))
    esc_v_perigee = t.esc_v(t.norm(r0_perigee))

    # Calculate Current Velocities
    v0_apogee_norm = t.norm(v0_apogee)
    v0_perigee_norm = t.norm(v0_perigee)

    # Calculate Escape Trajectory Vectors
    v0_apogee_escape = v0_apogee_normed * esc_v_apogee
    v0_perigee_escape = v0_perigee_normed * esc_v_perigee

    # Print Velocities and Differences
    print('Current Velocity at Apogee: \t%.2f km/s' % v0_apogee_norm)
    print('Circular Velocity at Apogee: \t%.2f km/s' % v_circ_apogee)
    print('Escape Velocity at Apogee: \t%.2f km/s' % esc_v_apogee)
    print('Delta V: \t%.2f km/s' % (esc_v_apogee - v0_apogee_norm))
    print('-------------------------------------------------')
    print('Current Velocity at Perigee: \t%.2f km/s' % v0_perigee_norm)
    print('Circular Velocity at Perigee: \t%.2f km/s' % v_circ_perigee)
    print('Escape Velocity at Perigee: \t%.2f km/s' % esc_v_perigee)
    print('Delta V: \t%.2f km/s' % (esc_v_perigee - v0_perigee_norm))

    # Create Perturbations Dictionary
    perts = null_perts()
    perts['N_Bodies'] = [pd.Moon]

    # Initial Conditions
    state0_original = r0_apogee.tolist() + v0_apogee.tolist()
    state0_apogee = r0_apogee.tolist() + v0_apogee_escape.tolist()
    state0_perigee = r0_perigee.tolist() + v0_perigee_escape.tolist()

    # Create Instances of Orbit Propagator
    op_original = OP(state0_original, tspan, dt, perts = perts, date0 = date0)
    op_apogee = OP(state0_apogee, tspan, dt, perts = perts, date0 = date0)
    op_perigee = OP(state0_perigee, tspan, dt, perts = perts, date0 = date0)
    op_rotated = OP(coes0_rotated, tspan, dt, coes = True, deg = True, perts = perts, date0 = date0)

    # Calculate and Plot COEs of Escape Trajectory from Perigee
    op_perigee.calculate_coes()
    op_perigee.plot_coes(rel = False, days = True)

    # Extract Moon State Vectors
    moon_traj = op_apogee.perts['N_Bodies'][0]['states'][:,:3]

    # Plot All
    t.plot_n_orbits([op_original.rs, op_rotated.rs, moon_traj, op_apogee.rs, op_perigee.rs],
                    ['Original Orbit', 'Rotated Orbit', 'Moon', 'Apogee Escape', 'Perigee Escape'],
                    show_plot = True)