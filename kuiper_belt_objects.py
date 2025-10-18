# Kuiper Belt Objects

import sys 
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import tools as t
import matplotlib.pyplot as plt
import spiceypy as spice
import spice_tools as st
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP

dt = 1000000
cb = pd.Sun

# Initial Date
date0 = '2020-01-01'

# Reference Frame
FRAME = 'ECLIPJ2000'

# Centre of Reference Frame
OBSERVER = '0'

# Header of Output CSV File
HEADER = 't, rx, ry, rz'

IDS = [999, 2010199, 2136108, 2136472]
NAMES = ['Neptune', 'Pluto', 'Chariklo', 'Haumea', 'Makemake']

if __name__ == '__main__':

    # Load SPICE Files
    spice.furnsh('spice_data/latest_leapseconds.tls')
    spice.furnsh('spice_data/de432s.bsp')
    spice.furnsh('spice_data/nh_plu017.bsp')
    spice.furnsh('spice_data/kbo_centaur_20131129.bsp')

    # Convert Start Date to Seconds past J2000
    et0 = spice.utc2et(date0)

    # Get Initial State of Neptune
    state0 = st.get_ephemeris_data('8', et0, FRAME, OBSERVER)

    # Propagate 1 Orbit of Neptune and Place in Array
    ephemerides = [OP(state0, '1', dt, cb = cb).rs]

    # Calculate Kuiper Belt Object Ephemerides
    for id_ in IDS:
        id_ = str(id_).strip()  # Clean the ID
        if id_:  # Skip Empty Strings
            ephemerides.append(OP(st.get_ephemeris_data(id_, et0, FRAME, OBSERVER), '1', dt, cb = cb).rs)

    if True:
        t.plot_n_orbits(ephemerides, labels = NAMES, cb = cb, title = 'Kuiper Belt Objects', el = 20, az = -45, AU = True, axes = 30, show_plot = True)

        