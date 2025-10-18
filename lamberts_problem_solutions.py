# Universal Variables Solution to Lambert's Problem

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
import lambert_tools as lt

dt = 1000

# Central Body
cb = pd.Sun

# Initial and Final Dates
date0 = '2005-12-01'
datef = '2006-03-01'

# Reference Frame
FRAME = 'ECLIPJ2000'

# Centre of Reference Frame
OBSERVER = 'SUN'

# Header of Output CSV File
HEADER = 't, rx, ry, rz'

if __name__ == '__main__':

    # Load SPICE Kernels
    spice.furnsh('spice_data/latest_leapseconds.tls')
    spice.furnsh('spice_data/de432s.bsp')

    # Convert Start and End Dates to Seconds Past J2000
    et0 = spice.utc2et(date0)
    etf = spice.utc2et(datef)

    # Calculate Transfer Time
    transfer_time = etf - et0

    # Time Arrays for Earth and Venus
    time_arr = np.linspace(et0, etf, 10000)

    # Calculate Eath and Venus State Vectors at Initial Time
    states_earth = st.get_ephemeris_data('EARTH', time_arr, FRAME, OBSERVER)
    states_venus = st.get_ephemeris_data('VENUS', time_arr, FRAME, OBSERVER)

    # Spacef=craft Initial Position Vector
    r0 = states_earth[0,:3]

    # Spacecraft Final Position Vector
    rf = states_venus[-1,:3]

    # Calculate Spacecraft Velocity Vectors via Lambert's Solution
    v0, vf = lt.lamberts_universal_variables(r0, rf, transfer_time, mu = cb['mu'])

    # Initial State Vector for Spacecraft
    state0_sc = r0.tolist() + v0.tolist()

    # Propagate Spacecraft Orbit
    op_sc = OP(state0_sc, transfer_time, dt, cb = cb)

    # t.write_csv(HEADER, [op_sc.ts, op_sc.rs], 'v_1 /earth_to_venus/spacecraft_traj.csv)
    # t.write_csv(HEADER, [time_arr, states_earth], 'v_1 /earth_to_venus/earth_traj.csv)
    # t.write_csv(HEADER, [time_arr, states_venus], 'v_1 /earth_to_venus/venus_traj.csv)

    # Plot
    t.plot_n_orbits([states_earth[:,:3], states_venus[:,:3], op_sc.rs], labels = ['Earth', 'Venus', 'Spacecraft'], cb = cb, show_plot = True, save_plot = False, title = 'Earth to Venus Transfer')


