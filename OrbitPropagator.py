# Orbit Propagator Script
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from scipy.integrate import ode
import planetary_data as pd
import tools as t
import spiceypy as spice
import spice_tools as st
import spice_data as sd
from time import time
from multiprocessing import Pool


style.use('dark_background')

hours = 3600.0
days = hours * 24



def null_perts():
    return {
            'J2': False,
            'J3': False,
            'J4': False,
            'J5': False,
            'J6': False,
            'J7': False,
            'Aero': False,
            'Cd': 0,
            'A': 0,
            'N_Bodies': [],
            'SRP': False,
            'A_SRP': 0,
            'SRP_Custom_Func': False,
            'CR': 0,
            'B': 0,
            'Oblateness': False,
            'Relativity': False,
            'Moon_Gravity': False,
            'Solar_Gravity': False,
            'Thrust': 0,
            'Thrust_Direction': 0,
            'Custom_Thrust_Function': False,
            'ISP': 0,
            'Rho': 0,
            'C20': 0,
            'Custom_Pert': False
        }

class OrbitPropagator:
    def __init__(self, state0, tspan, dt, coes = False, deg = True, mass0 = 0, date0 = 0, frame = 'J2000', cb = pd.Earth, perts = null_perts(), propagator = 'lsoda', sc = {}):
        # Check if Initial State is COE or RV
        if coes:
            self.r0, self.v0 = t.coes2rv(state0, deg = deg, mu = cb['mu'])
        else:
            self.r0 = np.array(state0[0:3])
            self.v0 = np.array(state0[3:])


        # Internally Store Parameters
        self.cb = cb
        self.propagator = propagator
        self.perts = perts
        self.tspan = tspan
        self.dt = dt
        self.mass0 = mass0
        self.date0 = date0
        self.frame = frame
        
        # Check if Time Span is One Period of n seconds
        if type(tspan) == str:
            self.tspan = float(tspan) * t.rv2period(self.r0, self.v0, self.cb['mu'])
        else:
            self.tspan = tspan
        


        # Miscellaneous Parameters
        self.n_steps = int(np.ceil(self.tspan / self.dt)) + 1
        self.ts = np.zeros((self.n_steps + 1, 1))
        self.y = np.zeros((self.n_steps + 1, 7))
        self.alts = np.zeros((self.n_steps + 1))
        self.step = 0
        self.spice_files_loaded = []

        # Initial Conditions
        self.y[0,:] = self.r0.tolist() + self.v0.tolist() + [self.mass0]
        self.alts[0] = t.norm(self.r0) - self.cb['radius']

        # Initialize Solver
        self.solver = ode(self.diff_eq)
        self.solver.set_integrator(self.propagator)
        self.solver.set_initial_value(self.y[0,:], 0)

        # Store Stop Conditions Dictionary
        self.stop_conditions_dict = sc

        # Define Dictionary to Map Internal Methods
        self.stop_conditions_map = {'max_alt': self.check_max_alt, 'min_alt': self.check_min_alt}

        # Create Stop Conditions Function List
        self.stop_condition_functions = [self.check_deorbit]

        # Fill in the Rest of the Stop Conditions
        for key in self.stop_conditions_dict:
            self.stop_condition_functions.append(self.stop_conditions_map[key])

        # Check if Loading in SPICE Data
        if self.perts['N_Bodies'] or self.perts['SRP']:

            # Load LeapSeconds Kernel
            spice.furnsh(sd.leap_seconds_kernel)

            # Add to List of Loaded SPICE Files
            self.spice_files_loaded.append(sd.leap_seconds_kernel)

            # Convert Start Date to Seconds after J2000
            self.start_time = spice.utc2et(self.date0)

            # Create Timespan Array in Seconds after J2000
            self.spice_tspan = np.linspace(self.start_time, self.start_time + self.tspan, self.n_steps)

            # If SRP, Get States of Sun
            if self.perts['SRP']:

                # Load SPICE Files for Given Body
                spice.furnsh(self.cb['spice_file'])

                # Add to List of Loaded SPICE Files
                self.spice_files_loaded.append(self.cb['spice_file'])

                # Calculate Central Body's Stats Throughout Enire Propagation wrt Sun
                self.cb['states'] = st.get_ephemeris_data(self.cb['name'], self.spice_tspan, self.frame, 'SUN')

                
        # Load Kernels for Each Body
        for body in self.perts['N_Bodies']:

            # If SPICE File Not Already Loaded, Load it
            if body['spice_file'] not in self.spice_files_loaded:

                # Load SPICE File
                spice.furnsh(body['spice_file'])

                # Append to List of Loaded SPICE Files
                self.spice_files_loaded.append(body['spice_file'])

            # Calculate Body States wrt Central Body
            body['states'] = st.get_ephemeris_data(body['name'], self.spice_tspan, self.frame, self.cb['name']) 

        # Check for Custom Thrust Function
        if self.perts['Custom_Thrust_Function']:
            self.thrust_func = self.perts['Custom_Thrust_Function']
        else:
            if self.perts['Thrust']:
                self.thrust_func = self.default_thrust_func

        
        self.propagate_orbit()

    # Check if Spacecraft has Deorbited
    def check_deorbit(self):
        if self.alts[self.step] < self.cb['deorbit_altitude']:
            print('Spacecraft deorbited after %.1f seconds' % self.ts[self.step])
            return False
        else:
            return True

    # Check if Maximum Altitude Exceeded
    def check_max_alt(self):
        if self.alts[self.step] > self.stop_conditions_dict['max_alt']:
            print('Spacecraft reached maximum altitude after %.1f seconds' % self.ts[self.step])
            return False
        else:
            return True

    # Check if Minimum Altitude Exceeded
    def check_min_alt(self):
        if self.alts[self.step] < self.stop_conditions_dict['min_alt']:
            print('Spacecraft reached minimum altitude after %.1f seconds' % self.ts[self.step])
            return False
        else:
            return True


    # Function Called at Each Time Step to Check Stop Conditions
    def check_stop_conditions(self):
        # For Each Stop Condition
        for sc in self.stop_condition_functions:
            # If returns False, Stop Propagation
            if not sc():
                return False
        else:
            return True

    
    # Propagate Orbit based on Given Conditions in init() Function
    def propagate_orbit(self):
        print('Propagating Orbit...')

        while self.solver.successful() and self.step < self.n_steps and self.check_stop_conditions():
            # Integrate Orbit           
            self.solver.integrate(self.solver.t + self.dt)
            self.step += 1
            # Extract Values from Solver Intance
            self.ts[self.step] = self.solver.t
            self.y[self.step] = self.solver.y

            # Calculate Altitude at this Time Step
            self.alts[self.step] = t.norm(self.solver.y[:3]) - self.cb['radius']


        # Extract Arrays at the Step Where Propagation Stopped
        self.ts = self.ts[0:self.step]
        self.rs = self.y[:self.step, :3]
        self.vs = self.y[:self.step ,3:6]
        self.masses = self.y[:self.step, 6]
        self.alts = (np.linalg.norm(self.rs, axis=1) - self.cb['radius']).reshape(self.step, 1)


    def diff_eq(self, _t, y):
        rx, ry, rz, vx, vy, vz, mass = y
        r = np.array([rx, ry, rz])
        v = np.array([vx, vy, vz])

        # Initialise dmdt
        dmdt = 0.0

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

        # Thrust Perturbations
        if self.perts['Thrust']:
            # Thrust Vector
            a += self.perts['Thrust Direction'] * t.normed(v) * self.perts['Thrust'] / mass / 1000.0

            # Derivative of Total Mass
            dmdt = -self.perts['Thrust'] / (self.perts['ISP'] * 9.81)


        # N-Body Perturbations
        for body in self.perts['N_Bodies']:

            # Get Vector Pointing from Central Body to Body
            r_cb2nb = body['states'][self.step, :3]

            # Get Vector Pointing from Satellite to Body
            r_sat2body = r_cb2nb - r

            # Nth Body Acceleration Vector
            a += body['mu'] * (r_sat2body / t.norm(r_sat2body)**3 - r_cb2nb / t.norm(r_cb2nb)**3)

        # Solar Radiation Pressure
        if self.perts['SRP']:

            # Vector Pointing from Sun to Spacecraft
            r_sun2sc = self.cb['states'][self.step, :3] + r

            # SRP Vector from Given Ephemeris Data
            a += (1 + self.perts['CR']) * pd.Sun['G1'] * self.perts['A_SRP'] / mass / t.norm(r_sun2sc)**3 * r_sun2sc
    
        return [vx, vy, vz, a[0], a[1], a[2], dmdt]


    # Custom Thrust Function (Velocity or Anti-Velocity Direction)
    def default_thrust_func(self, t_, state, _):
        return self.perts['Thrust Direction'] * t.normed(state[3:6]) * self.perts['Thrust'] / state[6] / 1000.0 # km/s^2


    # Calculate Classical Orbital Elements from Position and Velocity Vectors
    def calculate_coes(self, degrees = True, parallel = False):
        print('Calculating Classical Orbital Elements...')

        if parallel:
            start = time()
            states = t.parallel_states(self.y)
            pool = Pool()
            self.coes = np.array(pool.starmap(t.rv2coes, states))
            self.coes_rel = self.coes - self.coes[0,:]
            print(time() - start)

        else:
            start = time()

            # Preallocate Arrays
            self.coes = np.zeros((self.step, 6))
            self.coes_rel = np.zeros((self.step, 6))

            # Fill Arrays
            for n in range(self.step): 
                self.coes[n, :] = t.rv2coes(self.rs[n, :], self.vs[n, :], mu = self.cb['mu'], degrees = degrees, print_results = False)

            self.coes_rel = self.coes - self.coes[0,:]

            print(time() - start)


    # Plot Classical Orbital Elements vs Time
    def plot_coes(self, hours = False, days = False, show_plot = False, save_plot = False, title = 'Classical Orbital Elements', rel = True, figsize = (16, 8)):
        
        print('Plotting Classical Orbital Elements...')

        # Create Figure and Axes Instances
        fig, axs = plt.subplots(3, 2, figsize = figsize)
        # Figure Title
        fig.suptitle(title, fontsize = 16)
        # X Axis
        if hours:
            ts = self.ts / 3600.0
            xlabel = 'Time (hours)'
        elif days:
            ts = self.ts / 3600.0 / 24.0
            xlabel = 'Time (days)'
        else:
            ts = self.ts
            xlabel = 'Time (seconds)'


        # Check to Plot Relative COEs   
        if rel:
            coes = self.coes_rel
            title += ' (Relative)'
        else:
            coes = self.coes

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

        # Spread out Subplots
        plt.subplots_adjust(wspace = 0.3)

        if show_plot:
            plt.show()
        
        if save_plot:
            fig.savefig(title +'.png', dpi = 300)

    # Calculate Apoapse and Periapse
    def calculate_apoapse_periapse(self):
        # Define Empty Arrays
        self.apoapses = self.coes[:,0] * (1 + self.coes[:,1])
        self.periapses = self.coes[:,0] * (1 - self.coes[:,1])

    # Plot Apoapse and Periapse Over Time
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

    # 3D Plot of Orbit
    def plot_3d(self, show_plot = False, save_plot = False, title = None):
        fig = plt.figure(figsize=(16, 8))
        ax = fig.add_subplot(111, projection='3d')

        if title:
            plt.title(title)

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