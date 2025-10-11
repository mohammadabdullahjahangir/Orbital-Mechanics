# Using NASA Horizon Data to Track the Tesla Roadster in Space

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
 
HEADER = 't, rx, ry, rz'
DATE0 = '2018-02-08'
DATEF = '2025-10-10'

BODIES = ['MERCURY', 'VENUS', 'EARTH', 'MARS BARYCENTER', 'ROADSTER']
FRAME = 'ECLIPJ2000'
OBSERVER = 'SOLAR SYSTEM BARYCENTER'
STEPS = 1000

X = -1.112263266352074E+08; Y = 9.769962459386614E+07; Z = -9.295017561224103E+04
VX = -2.265237130832758E+01; VY =-2.523718038768873E+01; VZ = -7.067260463458140E-01
state0 = np.array([X, Y, Z, VX, VY, VZ])

if __name__ == '__main__':

    # Load Kernels for Solar System Ephemeris
    spice.furnsh('spice_data/latest_leapseconds.tls')
    spice.furnsh('spice_data/de432s.bsp')

    # Seconds Since J2000 for Initial Time
    t0 = spice.utc2et(DATE0)

    # Final Time in Seconds
    tf = spice.utc2et(DATEF)

    # tspan for Orbit Propagator
    tspan = tf - t0

    # Create Time Array for Ephemeris Data for all Bodies (Assuming all are the same)
    times = st.tc2array([t0,tf], STEPS)

    # Create Empty List for all Ephemeris Data
    rs = []

    # For each Body in Solar System
    for name in BODIES[:-1]:
        
        # Add Ephemeris Data to List
        rs.append(st.get_ephemeris_data(name, times, FRAME, OBSERVER))

        t.writecsv(HEADER, t.joinarrs([times, rs[-1]]), 'roadster/' + name + '_traj.csv')

    # Create Instance of Orbit Propagator for Tesla Roadster
    op0 = OP(state0, tspan, 1000, propagator = 'lsoda', cb = pd.Sun)

    # Write Trajectroy to CSV 
    op0.write_traj('roadster/roadster_traj.csv')

    # Add Position Vectors to List
    rs.append(op0.rs)

    # Plot 
    t.plot_n_orbits(rs, labels = BODIES, cb = pd.Sun, save_plot = True)

    