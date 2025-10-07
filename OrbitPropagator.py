import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from scipy.integrate import ode
import planetary_data as pd
import tools as t
#style.use('dark_background')

class OrbitPropagator:
    def __init__(self, state0, r0, v0, tspan, dt, coes = False, deg = True, cb = pd.Earth):
        if coes:
            self.r0, self.v0 = t.coes2rv(state0, deg = deg, mu = cb['mu'])
        else:
            self.r0 = state0[0:3]
            self.v0 = state0[3:]
        self.cb = cb
        self.tspan = tspan
        self.dt = dt
        
        self.n_steps = int(np.ceil(self.tspan / self.dt))
        self.ys = np.zeros((self.n_steps, 6))
        self.ts = np.zeros((self.n_steps, 1))

        self.y0 = self.r0.tolist() + self.v0.tolist()
        self.ys[0] = self.y0
        self.step = 1

        self.solver = ode(self.two_body)
        self.solver.set_integrator('dopri5')
        self.solver.set_initial_value(self.y0, 0)
        self.solver.set_f_params(self.cb['mu'])

        self.propagate_orbit()

    def propagate_orbit(self):
        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t + self.dt)
            self.ts[self.step] = self.solver.t
            self.ys[self.step] = self.solver.y
            self.step += 1

        self.rs = self.ys[:, 0:3]
        self.vs = self.ys[:,3:]
        
        
    def two_body(self, t, y, mu):
        rx, ry, rz, vx, vy, vz = y
        r = np.array([rx, ry, rz])
        v = np.array([vx, vy, vz])
        r_norm = np.linalg.norm(r)

        ax, ay, az = -r * self.cb['mu'] / r_norm**3
        return [vx, vy, vz, ax, ay, az]


    def plot3d(self):
        ax = plt.axes(projection='3d')
        ax.plot(self.rs[:, 0], self.rs[:, 1], self.rs[:, 2], 'w', label ='Trajectory')
        ax.plot(self.rs[0, 0], self.rs[0, 1], self.rs[0, 2], 'wo', markersize=10, label ='Starting Position')

        r_plot = self.cb['radius']

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

        max_value = np.max(np.abs(self.rs))

        ax.set_xlim([-max_value, max_value])
        ax.set_ylim([-max_value, max_value])
        ax.set_zlim([-max_value, max_value])
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
        ax.set_aspect('equal')
        plt.legend(['Trajectory', 'Starting Position'])

