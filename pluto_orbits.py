# Pluto's Planetary Status

import sys 
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import tools as t
import spiceypy as spice
import spice_tools as st

# Initial Date
date0 = '2020-12-01'

# Reference Frame
FRAME = 'ECLIPJ2000'

# Centre of Referenece Frame
OBSERVER = 'PLUTO BARYCENTER'

# Header of Output CSV File
HEADER = 't, rx, ry, rz, vx, vy, vz'

NAMES = ['PLUTO', 'CHARON', 'NIX', 'HYDRA']

if __name__ == '__main__':

    for name in NAMES:
        os.makedirs(f'pluto/{name}', exist_ok=True)

    # Load SPICE Kernels
    spice.furnsh('spice_data/latest_leapseconds.tls')
    spice.furnsh('spice_data/de432s.bsp')
    spice.furnsh('spice_data/nh_plu017.bsp')

    # Convert Start Date to Seconds Past J2000
    et0 = spice.utc2et(date0)

    # End Time after ~2 Orbits of Charon
    etf = et0 + 1808650

    # Get Linspace of Times
    times = st.tc2array([et0, etf], 1000)

    # Caclulate Pluto System Bodies Ephemerides
    ephemerides = []

    for name in NAMES:

        # Get Initial State Vector
        ephemerides.append(st.get_ephemeris_data(name, times, FRAME, OBSERVER))

        # Write Output CSV
        eph = ephemerides[-1]
        data_to_write = np.column_stack([times, eph])
        t.write_csv(HEADER, data_to_write, filename = 'pluto/' + name + '/traj.csv')

    
    t.plot_n_orbits(ephemerides, labels = NAMES, title = 'Pluto System', axes = 40000, el = 30, az = -30, show_plot = True, dpi = 300)
