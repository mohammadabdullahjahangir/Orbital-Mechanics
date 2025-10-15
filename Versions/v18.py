# Low Thrust Escape Velocities

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
tspan = 3600 * 24 * 20 # 20 days
dt = 100

tspan1 = 3600 * 24 * 150

date0 = '2020-03-28'

if __name__ == '__main__':

    # Define Perturbations Dictionary for Low Thrust Model
    perts = null_perts()
    perts['Thrust'] = 0.327
    perts['ISP'] = 4300
    perts['Thrust_Direction'] = 1
    perts['N_Bodies'] = [pd.Moon]

    # Define Stop Conditions for Low Thrust Trajectory
    sc = {'escape_velocity': True}

    # Initial Mass of Spacecraft
    mass0 = 50

    # Calculate Initial State Vector
    state0 = [cb['radius'] + 800, 0.01, 5.0, 0.0, 0.0, 0.0]

    # Create Instance of Orbit Propagator
    op_low_thrust = OP(state0, tspan, dt, deg = True, coes = True, mass0 = mass0, perts = perts, sc = sc, date0 = date0)

    # Plot Altitudes over Time
    op_low_thrust.plot_alts(hours = True)

    # Plot Trajectory
    op_low_thrust.plot_3d(show_plot = True)

    # Calculate and Plot COEs
    op_low_thrust.calculate_coes()
    op_low_thrust.plot_coes(hours = True, rel = False)

    # Calculate and Plot Apoapse and Periapse
    op_low_thrust.calculate_apoapse_periapse()
    op_low_thrust.plot_apoapse_periapse(show_plot = True, hours = True)

    # Create new Instance Perturbations Dictionary 
    perts = null_perts()
    perts['N_Bodies'] = [pd.Moon]

    # Create New Instances of Orbit Propagator with Initial Conditions as Ending of Previous
    op1 = OP(op_low_thrust.y[-1,:-1], tspan1, dt, perts = perts)
    op2 = OP(op_low_thrust.y[-100,:-1], tspan1, dt, perts = perts)

    # Extract Moon State Vectors
    moon_traj = op1.perts['N_Bodies'][0]['states'][:,:3]

    print("op1.y shape:", op1.y.shape)
    print("op1.rs shape:", op1.rs.shape)
    print("op2.y shape:", op2.y.shape)
    print("op2.rs shape:", op2.rs.shape)

    # Plot All
    t.plot_n_orbits([op_low_thrust.rs, op1.rs, op2.rs, moon_traj], ['Low Thrust', 'Escape Coast', 'Non-Escape Coast', 'Moon'])

