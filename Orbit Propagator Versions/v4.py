# Orbit propagation Test Script
import sys 
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import style
from scipy.integrate import ode
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP

style.use('dark_background')

cb = pd.Earth

if __name__ == '__main__':
    r_mag = cb['radius'] + 500.0  # Initial orbit radius (km)
    v_mag = np.sqrt(cb['mu'] / r_mag)    # Initial orbit velocity (km/s)
    r0 = [r_mag, 0.0, 0.0]  # Initial position vector (km)
    v0 = [0.0, v_mag, 0.0]  # Initial velocity vector (km/s)

    tspan = 24 * 60 * 60  # Total time span (s)
    dt = 100.0

    n_steps = int(np.ceil(tspan / dt))

    ys = np.zeros((n_steps, 6))
    ts = np.zeros((n_steps, 1))

    y0 = r0 + v0
    ys[0] = np.array(y0)
    step = 1

    op = OP(r0, v0, tspan, dt)
    op.propagate_orbit()
    op.plot3d()