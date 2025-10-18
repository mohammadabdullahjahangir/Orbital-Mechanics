# Solar System Orbits with NASA SPICE Files

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

# Total Steps of Epehemeris Data
STEPS = 100000

# Reference Frame 
FRAME = 'ECLIPJ2000'

# Observer
OBSERVER = 'SUN'

if __name__ == '__main__':

    # Load SPICE Kernels
    
    spice.furnsh('spice_data/latest_leapseconds.tls')
    spice.furnsh('spice_data/de432s.bsp')

    # Get Object IDs, Names, and Time Coverages in SPK File
    ids, names, tcs_sec, tcs_cal = st.get_objects('spice_data/de432s.bsp', display=True)

    # Only Include BARYCENTERS
    names = [f for f in names if 'BARYCENTER' in f]

    # Create Time Array for Epheremis Data
    times = st.tc2array(tcs_sec[0], STEPS)

    # Create Empty List for all Ephemeris Data
    rs = []

    # For each Body in Solar System
    for name in names:
        # Add Ephemeris Data to List
        rs.append(st.get_ephemeris_data(name, times, FRAME, OBSERVER))

    t.plot_n_orbits(rs, names, show_plot = True, AU = True, cb = pd.Sun, figsize = (16,16))


    
