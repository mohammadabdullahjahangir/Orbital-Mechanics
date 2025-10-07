import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from scipy.integrate import ode
import planetary_data as pd
import tools as t
style.use('dark_background')


def null_perts():
    return {
            'J2': False,
            'Aero': False,
            'Oblateness': False,
            'Moon_Gravity': False,
            'Solar_Gravity': False,

        }

class OrbitPropagator:
    def __init__(self, state0, r0, v0, tspan, dt, coes = False, deg = True, cb = pd.Earth, perts = null_perts()):
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

        self.solver = ode(self.diff_eq)
        self.solver.set_integrator('dopri5')
        self.solver.set_initial_value(self.y0, 0)
        self.solver.set_f_params(self.cb['mu'])
        
        # Define Perturbations Dictionary
        self.perts = perts
        self.mass = perts.get('mass', 1000.0)  # default mass = 1000 kg


        self.propagate_orbit()

    def propagate_orbit(self):
        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t + self.dt)
            self.ts[self.step] = self.solver.t
            self.ys[self.step] = self.solver.y
            self.step += 1

        # Extract Arrays at the Step Where Integration Ended
        self.ts = self.ts[0:self.step]
        self.rs = self.ys[:self.step, :3]
        self.vs = self.ys[:self.step ,3:]
        self.alts = (np.linalg.norm(self.rs, axis=1) - self.cb['radius']).reshape(self.step, 1)

    def diff_eq(self, _t, y, mu):
        rx, ry, rz, vx, vy, vz = y
        r = np.array([rx, ry, rz])
        v = np.array([vx, vy, vz])

        # Magnitude of Position Vector
        r_norm = np.linalg.norm(r)

        # Two Body Acceleration
        a = -r*self.cb['mu']/r_norm**3

        # Oblateness Perturbations
        if self.perts['Oblateness']:
            z2 = r[2]**2
            r2 = r_norm**2
            tx = r[0] / r_norm * (5 * z2 / r2 - 1)
            ty = r[1] / r_norm * (5 * z2 / r2 - 1)
            tz = r[2] / r_norm * (5 * z2 / r2 - 3)
            a += 1.5 * self.cb['J2'] * self.cb['mu'] * self.cb['radius']**2 / r_norm**4 * np.array([tx, ty, tz])

        # Aerodynamic Drag Perturbation
        if self.perts['Aero']:
            # Calculate Altitude and Air Density
            z = r_norm - self.cb['radius']
            rho = t.calc_atmospheric_density(z)

            # Calculate Motion of Spacecraft Relative to Atmosphere
            v_rel = v - np.cross(self.cb['atm_rot_vector'], r)
            drag = -v_rel * 0.5 * rho * np.linalg.norm(v_rel) * self.perts['Cd'] * self.perts['A'] / self.mass
            a += drag
        
        # J2 Perturbation
        if self.perts['J2']:
            z2 = r[2]**2
            r2 = r_norm**2
            tx = r[0]/r_norm*(5*z2/r2 - 1)
            ty = r[1]/r_norm*(5*z2/r2 - 1)
            tz = r[2]/r_norm*(5*z2/r2 - 3)

            a_j2 = 1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/r_norm**4 * np.array([tx, ty, tz])
            a += a_j2

        
        return [vx, vy, vz, a[0], a[1], a[2]]


    def calculate_coes(self, degrees = True):
        print('Calculating Classical Orbital Elements...')
        self.coes = np.zeros((self.n_steps, 6))
        for n in range(self.n_steps): 
            self.coes[n, :] = t.rv2coes(self.rs[n, :], self.vs[n, :], mu = self.cb['mu'], degrees = degrees, print_results = False)

    def plot_coes(self, hours = False, show_plot = False, save_plot = False, title = 'Classical Orbital Elements', figsize = (16, 8)):
        
        print('Plotting Classical Orbital Elements...')

        # Create Figure and Axes Instances
        fig, axs = plt.subplots(3, 2, figsize = figsize)
        # Figure Title
        fig.suptitle(title, fontsize = 16)
        # X Axis
        if hours:
            ts = self.ts / 3600.0
            xlabel = 'Time (hours)'
        elif self.days:
            ts = self.ts / 3600.0 / 24.0
            xlabel = 'Time (days)'
        else:
            ts = self.ts
            xlabel = 'Time (seconds)'

        # Plot True Anomaly
        axs[0,0].plot(self.ts, self.coes[:,3])
        axs[0,0].set_title('True Anomaly vs Time')
        axs[0,0].grid(True)
        axs[0,0].set_xlabel(xlabel)
        axs[0,0].set_ylabel('True Anomaly (degrees)')
        
        # Plot Semi-Major Axis
        axs[1,0].plot(self.ts, self.coes[:,0])
        axs[1,0].set_title('Semi-Major Axis vs Time')
        axs[1,0].grid(True)
        axs[1,0].set_xlabel(xlabel)
        axs[1,0].set_ylabel('Semi-Major Axis (km)')

        # Plot Eccentricity
        axs[0,1].plot(self.ts, self.coes[:,1])
        axs[0,1].set_title('Eccentricity vs Time')
        axs[0,1].grid(True)
        axs[0,1].set_xlabel(xlabel)
        axs[0,1].set_ylabel('Eccentricity')

        # Plot Argument of Periapse
        axs[2,0].plot(self.ts, self.coes[:,4])
        axs[2,0].set_title('Argument of Periapse vs Time')
        axs[2,0].grid(True)

        # Plot Inclination
        axs[1,1].plot(self.ts, self.coes[:,2])
        axs[1,1].set_title('Inclination vs Time')
        axs[1,1].grid(True)
        axs[1,1].set_xlabel(xlabel)
        axs[1,1].set_ylabel('Inclination (degrees)')

        # Plot Right Ascension of the Ascending Node
        axs[2,1].plot(self.ts, self.coes[:,5])
        axs[2,1].set_title('RAAN vs Time')
        axs[2,1].grid(True)

        plt.tight_layout(rect=[0, 0, 1, 0.97])

        if show_plot:
            plt.show()
        
        if save_plot:
            fig.savefig(title +'.png', dpi = 300)

    
    def calculate_apoapse_periapse(self):
        # Define Empty Arrays
        self.apoapses = self.coes[:,0] * (1 + self.coes[:,1])
        self.periapses = self.coes[:,0] * (1 - self.coes[:,1])

    def plot_apoapse_periapse(self, hours = False, days = False, show_plot = False, save_plot = False, title = 'Apoapse and Periapse', figsize = (20, 10)):
        fig = plt.figure(figsize = figsize)

        if hours:
            ts = self.ts / 3600.0
            x_unit = 'Hours'
        elif days:
            ts = self.ts / 3600.0 / 24.0
            x_unit = 'Days'
        else:
            ts = self.ts
            x_unit = 'Seconds'

        plt.title(title)
        plt.xlabel('Time (%s)' % x_unit)
        plt.ylabel('Altitude (km)')
        plt.plot(ts, self.apoapses, 'b', label = 'Apoapse')
        plt.plot(ts, self.periapses, 'r', label = 'Periapse')
        plt.legend()
        plt.grid(True)

        if show_plot:
            plt.show()
        
        if save_plot:
            fig.savefig(title +'.png', dpi = 300)
        

        

    # Plot Altitude Over Time
    def plot_alts(self, hours = False, days = False, show_plot = False, save_plot = False, title = 'Radial Distance vs Time', figsize = (20, 10), dpi = 500):
        if hours:
            ts = self.ts / 3600.0
            x_unit = 'Hours'
        elif days:
            ts = self.ts / 3600.0 / 24.0
            x_unit = 'Days'
        else:
            ts = self.ts
            x_unit = 'Seconds'

        plt.figure(figsize = figsize)
        plt.plot(ts, self.alts, 'w')
        plt.grid(True)
        plt.xlabel('Time (%s)' % x_unit)
        plt.ylabel('Altitude (km)')
        plt.title(title)
        plt.show()

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title + '.png', dpi = dpi)


    def plot_3d(self, show_plot = False, save_plot = False):
        fig = plt.figure(figsize=(16, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.rs[:, 0], self.rs[:, 1], self.rs[:, 2], 'w', label ='Trajectory')
        ax.plot(self.rs[0, 0], self.rs[0, 1], self.rs[0, 2], 'wo', markersize = 10, label ='Starting Position')

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
        plt.show()