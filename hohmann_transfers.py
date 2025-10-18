# Hohmann Transfers

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
from Spacecraft import Spacecraft as SC

cb = pd.Earth

if __name__ == '__main__':

    # Initial and Final Altitudes
    r0 = 1000
    r1 = 10000

    # COEs of Initial and Final Orbits
    coes0 = [cb['radius'] + r0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coes1 = [cb['radius'] + r1, 0.0, 0.0, 180.0, 0.0, 0.0]

    # Rotated Orbital Elements
    coes0_rot = [cb['radius'] + r0, 0.0, 30.0, 0.0, 100.0, 20.0]
    coes1_rot = [cb['radius'] + r1, 0.0, 30.0, 180.0, 100.0, 20.0]

    # Call Scalar Hohmann Transfer Function
    delta_vs_scalar, t_transfer_scalar = t.hohmann_transfer_scalars(r0, r1)

    # Print to Terminal
    print('Delta V 0: \t %.3f (km/s)' % delta_vs_scalar[0])
    print('Delta V 1: \t %.3f (km/s)' % delta_vs_scalar[1])

    sc0, sc1, sc_transfer, delta_vs = t.hohmann_transfer(coes0 = coes0, coes1 = coes1, propagate = True)

    sc0_rot, sc1_rot, sc_transfer_rot, delta_vs_rot = t.hohmann_transfer(coes0 = coes0_rot, coes1 = coes1_rot, propagate = True)

    t.plot_n_orbits([sc0.rs, sc1.rs, sc_transfer.rs, sc0_rot.rs, sc1_rot.rs, sc_transfer.rs], labels = ['Initial', 'Final', 'Transfer', 'Initial_', 'Final_', 'Transfer_'], show_plot = True, az = -45, el = 0.0, axes = 9500, no_axes = True, title = '')


