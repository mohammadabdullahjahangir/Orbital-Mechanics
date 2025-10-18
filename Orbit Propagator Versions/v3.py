# Two Body Problem Propagator
import sys 
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode
from matplotlib import style

style.use('dark_background')


# Earth's gravitational parameter (km^3/s^2)
earth_radius = 6371.0
earth_mu = 398600.4418  

def plot(r):
    ax = plt.axes(projection='3d')
    ax.plot(r[:, 0], r[:, 1], r[:, 2], 'k')
    ax.plot([0], [0], [0], 'ko', markersize=10)

    r_plot = earth_radius

    # Create a sphere for Earth
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = r_plot * np.cos(_u) * np.sin(_v)
    _y = r_plot * np.sin(_u) * np.sin(_v)
    _z = r_plot * np.cos(_v)
    ax.plot_surface(_x, _y, _z, color='b', alpha=0.3)

    l = r_plot * 2.0
    x, y, z = [0, 0, 0], [0, 0, 0], [0, 0, 0]
    u, v, w = [1, 0, 0], [0, 1, 0], [0, 0, 1]
    ax.quiver(x, y, z, u, v, w, length=l, color='r', arrow_length_ratio=0.1)

    max_value = np.max(np.abs(r))

    ax.set_xlim([-max_value, max_value])
    ax.set_ylim([-max_value, max_value])
    ax.set_zlim([-max_value, max_value])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('equal')
    plt.legend(['Trajectory', 'Starting Position'])
    plt.show()


# 2-body equations of motion
def two_body(t, y, mu):
    rx, ry, rz, vx, vy, vz = y
    r = np.array([rx, ry, rz])
    r_norm = np.linalg.norm(r)

    ax, ay, az = -r * mu / r_norm**3
    return [vx, vy, vz, ax, ay, az]


if __name__ == '__main__':

    r_mag = earth_radius + 500.0  # Initial orbit radius (km)
    v_mag = np.sqrt(earth_mu / r_mag)    # Initial orbit velocity (km/s)
    r0 = [r_mag, 0.0, 0.0]  # Initial position vector (km)
    v0 = [0.0, v_mag, 0.0]  # Initial velocity vector (km/s)

    tspan = 100 * 60  # Total time span (s)
    dt = 100.0

    n_steps = int(np.ceil(tspan / dt))

    ys = np.zeros((n_steps, 6))
    ts = np.zeros((n_steps, 1))

    y0 = r0 + v0
    ys[0] = np.array(y0)
    step = 1

    solver = ode(two_body)
    solver.set_integrator('dopri5')
    solver.set_initial_value(y0, 0)
    solver.set_f_params(earth_mu)

    while solver.successful() and step < n_steps:
        solver.integrate(solver.t + dt)
        ts[step] = solver.t
        ys[step] = solver.y
        step += 1

    rs = ys[:, 0:3]

    plot(rs)    

