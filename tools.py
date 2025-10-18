# Tools Script
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import planetary_data as pd
import math as m
from datetime import datetime, timedelta
from Spacecraft import Spacecraft as SC

d2r = np.pi / 180.0  # Degrees to Radians
r2d = 180.0 / np.pi  # Radians to Degrees


km2AU = 1.4959787e8 # km to AU
figsize = (16, 8)


def norm(v):
    return np.linalg.norm(v)

def normed(v):
    return np.array(v) / norm(v)

# Write a header and data array to a CSV file
def write_csv(header, data, filename):
    """Write a header and data array to a CSV file."""
    # ensure numpy array
    data = np.array(data)
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)

# Join multiple 1D or 2D arrays column-wise, handling shape mismatches
def joinarrs(arrs):
    # Ensure all arrays are numpy arrays
    arrs = [np.atleast_2d(a).T if a.ndim == 1 else a for a in arrs]
    
    # Check that all arrays have the same number of rows
    n_rows = [a.shape[0] for a in arrs]
    if len(set(n_rows)) != 1:
        raise ValueError(f"Arrays have different lengths: {n_rows}")
    
    return np.hstack(arrs)

def plot_n_orbits(rs, labels, cb = pd.Earth, show_plot = False, save_plot = False, axes = False, AU = False, ER = False, figsize = (16,8), title = None, az = -60, el = 30, dpi = 300, output_dir = None, no_axes = False):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection='3d')

    # Set viewing angle
    ax.view_init(elev=el, azim=az)

    # Handle axes visibility
    if no_axes:
        ax.set_axis_off()

    if title is None:
        plt.title(title)
        
    # Plot Trajectories
    max_value = 0
    n = 0
    cs = ['c', 'orange', 'b', 'r', 'purple', 'brown', 'cyan', 'pink', 'gray']

    
    # Convert rs to list if needed and process
    rs_plot = []
    for r in rs:
        r_copy = r.copy()
        if AU:
            r_copy = r_copy / km2AU
        elif ER:
            r_copy = r_copy / pd.Earth['radius']
        rs_plot.append(r_copy)
    
    # Calculate max value from all orbits
    for r in rs_plot:
        max_ = np.max(np.abs(r))
        if max_ > max_value:
            max_value = max_
    
    # Plot each orbit
    n = 0
    for r in rs_plot:
        if labels is None:
            label0 = ''
        else:
            label0 = labels[n]
        
        # Plot Trajectory and Initial Position
        ax.plot(r[:, 0], r[:, 1], r[:, 2], c=cs[n % len(cs)], label=label0, linewidth=1)
        ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], 'o', c=cs[n % len(cs)], markersize=5)
        n += 1
    

    # Radius of Central Body
    r_plot = cb['radius']
    if AU:
        r_plot = r_plot / km2AU
    elif ER:
        r_plot = r_plot / pd.Earth['radius']
    
    # Plot Central Body (Sun)
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = r_plot * np.cos(_u) * np.sin(_v)
    _y = r_plot * np.sin(_u) * np.sin(_v)
    _z = r_plot * np.cos(_v)
    ax.plot_surface(_x, _y, _z, color='yellow', alpha=0.9)
    
    if axes:
        max_value = axes
    
    if AU:
        ax.set_xlabel('X (AU)')
        ax.set_ylabel('Y (AU)')
        ax.set_zlabel('Z (AU)')
    else:
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
    
    ax.set_xlim([-max_value, max_value])
    ax.set_ylim([-max_value, max_value])
    ax.set_zlim([-max_value, max_value])
    ax.set_box_aspect([1, 1, 1])
    
    # Set black background for better visibility
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')
    
    plt.legend(loc='upper left', fontsize=8)
    
    if show_plot:
        plt.show()
    if save_plot:
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            filepath = os.path.join(output_dir, f"{title or 'orbits'}.png")
        else:
            filepath = f"{title or 'orbits'}.png"
        fig.savefig(filepath, dpi=dpi)
        plt.savefig('n_orbits.png', dpi=300, facecolor='black')

# Define Hohmann Transfer Scalars
def hohmann_transfer_scalars(r0, r1, altitude = True, cb = pd.Earth):

    # if Passing in Altitude (Not Semi-Major Axis)
    if altitude:
        # Add Central Body Radius to r Values
        r0 += cb['radius']
        r1 += cb['radius']

    # Calculate Semi-Major Axis of Transfer Orbit
    a_transfer = (r0 + r1) / 2.0

    # Calculate Velocities of Ciruclar Orbits
    v_circ_init = m.sqrt(cb['mu'] / r0)
    v_circ_final = m.sqrt(cb['mu'] / r1)

    # Calculate Initial and Final Transfer Orbit Velocities (vis-viva equation)
    v0_transfer = m.sqrt(cb['mu'] * (2 / r0 - 1 / a_transfer))
    v1_transfer = m.sqrt(cb['mu'] * (2 / r1 - 1 / a_transfer))

    # Calculate Transfer Time
    t_transfer = m.pi * m.sqrt(a_transfer ** 3 / cb['mu'])

    # Calculate DeltaV Values
    delta_vs = [abs(v0_transfer - v_circ_init), abs(v_circ_final - v1_transfer)]

    return delta_vs, t_transfer

# Define Hohmann Transfer
def hohmann_transfer(r0 = 0, r1 = 0, cb = pd.Earth, coes0 = [], coes1 = [], altitude = True, propagate = False, dt = 100, output_dir = '', names = ['Initial', 'Final', 'Transfer'], write_output = False):
    
    # Check if COEs Passed In 
    if coes0:
        # Extract r0 and r1 Values
        r0 = coes0[0]
        r1 = coes1[0]

    # If Passing in Altitude (Not Semi-Major Axis)
    elif altitude:
        # Add Central Body Radius to r Values
        r0 += cb['radius']
        r1 += cb['radius']

    # Calculate Semi-Major Axis of Transfer Orbit
    a_transfer = (r0 + r1) / 2.0

    # Calculate Velocities of Ciruclar Orbits
    v_circ_init = m.sqrt(cb['mu'] / r0)
    v_circ_final = m.sqrt(cb['mu'] / r1)

    # Calculate Initial and Final Transfer Orbit Velocities (vis-viva equation)
    v0_transfer = m.sqrt(cb['mu'] * (2 / r0 - 1 / a_transfer))
    v1_transfer = m.sqrt(cb['mu'] * (2 / r1 - 1 / a_transfer))

    # Calculate Transfer Time
    t_transfer = m.pi * m.sqrt(a_transfer ** 3 / cb['mu'])

    # Calculate DeltaV Values
    delta_vs = [abs(v0_transfer - v_circ_init), abs(v_circ_final - v1_transfer)]

    # Propagate Orbits
    if propagate:
        # If COEs not Passed In
        if not coes0:

            # create COEs List
            coes0 = [r0, 0.0, 0.0, 0.0, 0.0, 0.0]
            coes1 = [r1, 0.0, 0.0, 180.0, 0.0, 0.0]

        # Calculate Eccentricity of Transfer Orbit
        e_transfer = 1 - r0 / a_transfer

        # COEs for Transfer Orbit
        coes_transfer = [a_transfer, e_transfer, coes0[2], 0.0, coes0[4], coes0[5]]

        # Calculate Periods of Initial and Final Orbits
        T0 = 2 * np.pi * (r0 ** 3 / cb['mu']) ** 0.5
        T1 = 2 * np.pi * (r1 ** 3 / cb['mu']) ** 0.5

        # Create Spacecraft Instances and Propagate Orbits
        sc0 = SC(coes0, T0, dt, coes = True, output_dir = os.path.join(output_dir, names[0]))
        sc1 = SC(coes1, T1, dt, coes = True, output_dir = os.path.join(output_dir, names[1]))
        sc_transfer = SC(coes_transfer, t_transfer, dt, coes = True, output_dir = os.path.join(output_dir, names[0]))

        # Write Output Files
        if write_output:
            sc0.write_traj()
            sc1.write_traj()
            sc_transfer.write_traj()
        
        return sc0, sc1, sc_transfer, delta_vs

    return delta_vs, t_transfer


# Define Escape Velocity
def esc_v(r, mu = pd.Earth['mu']):
    return m.sqrt( 2 * mu / r)

# Convert Classical Orbital Elements to Position and Velocity Vectors
def coes2rv(coes, deg=False, mu=pd.Earth['mu']):
    if deg:
        a, e, i, ta, aop, raan = coes
        i *= d2r
        ta *= d2r
        aop *= d2r
        raan *= d2r
    else:
        a, e, i, ta, aop, raan = coes
    
    E = ecc_anomaly([ta, e], 'tae')
    r_norm = a * (1 - e**2) / (1 + e * np.cos(ta))
    
    # Calculate Position and Velocity in Perifocal Frame
    r_perifocal = r_norm * np.array([m.cos(ta), m.sin(ta), 0.0])
    v_perifocal = m.sqrt(mu * a) / r_norm * np.array([-m.sin(E), m.sqrt(1 - e**2) * m.cos(E), 0.0])
    
    # Rotation Matrix from Perifocal to ECI Frame
    perif2eci = np.transpose(eci2perif(raan, aop, i))
    
    # Calculate Position and Velocity in ECI Frame
    r = np.dot(perif2eci, r_perifocal)
    v = np.dot(perif2eci, v_perifocal)
    
    return r, v

def rv2period(r, v, mu):
    # Calculate Specific Orbital Energy
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)
    
    # Specific Energy
    epsilon = v_mag**2 / 2 - mu / r_mag
    
    # Semi-Major Axis
    a = -mu / (2 * epsilon)
    
    # Orbital Period
    T = 2 * np.pi * np.sqrt(a**3 / mu)
    
    return T

def rv2coes(r, v, mu = pd.Earth['mu'], degrees = False, print_results = False):
    # Magnitudes of Position Vector
    r_norm = norm(r)

    # Specific Angular Momentum
    h = np.cross(r, v)
    h_norm = norm(h)

    # Inclination
    i = m.acos(h[2] / h_norm)

    # Eccentricity Vector
    e = (1 / mu) * ((norm(v)**2 - mu / r_norm) * r - np.dot(r, v) * v)

    # Eccentricity Scalar
    e_norm = norm(e)

    # Node Line
    N = np.cross([0, 0, 1], h)
    N_norm = norm(N)

    # Right Ascension of the Ascending Node
    raan = m.acos(N[0] / N_norm)
    if N[1] < 0:
        raan = 2 * np.pi - raan

    # Argument of Perigee
    aop = m.acos(np.dot(N, e) / (N_norm * e_norm))
    if e[2] < 0:
        aop = 2 * np.pi - aop

    # True Anomaly
    ta = m.acos(np.dot(e, r) / (e_norm * r_norm))
    if np.dot(r, v) < 0:
        ta = 2 * np.pi - ta
    
    # Semi-Major Axis
    a = r_norm*(1 + e_norm * m.cos(ta)) / (1 - e_norm**2)

    if print_results:
        print('a', a)
        print('e', e_norm)
        print('i', i * r2d)
        print('raan', raan * r2d)
        print('aop', aop * r2d)
        print('ta', ta * r2d)

    if degrees:
        return [a, e_norm, i * r2d, ta * r2d, aop * r2d, raan * r2d]
    else:
        return [a, e_norm, i, ta, aop, raan]

# Calculate Atmospheric Density at Altitude z (km)
def calc_atmospheric_density(z):
    rhos, zs = find_rho_z(z)
    if rhos[0] == 0:
        return 0.0
    
    Hi = -(zs[1] - zs[0]) / m.log(rhos[1] / rhos[0])

    return rhos[0] * m.exp(-(z - zs[0]) / Hi)

# Find Endpoints of Altitude and Density Surrounding Input Altitude
def find_rho_z(z, zs = pd.Earth['zs'], rhos = pd.Earth['rhos']):
    if not 1.0 < z < 1000.0:
        return [0.0, 0.0], [0.0, 0.0]

    # Find the Two Points Surrounding the Given Input Altitude
    for n in range(len(rhos) - 1):
        if zs[n] < z < zs[n + 1]:
            return [rhos[n], rhos[n + 1]], [zs[n], zs[n + 1]]
        
    # If Altitude is Outside Range of Table
    return [0.0, 0.0], [0.0, 0.0]


# Rotation Matrix from ECI to Perifocal Frame
def eci2perif(raan, aop, i):
    row0 = [-m.sin(raan) * m.cos(i) * m.sin(aop) + m.cos(raan) * m.cos(aop), 
            m.cos(raan) * m.cos(i) * m.sin(aop) + m.sin(raan) * m.cos(aop), 
            m.sin(i) * m.sin(aop)]
    row1 = [-m.sin(raan) * m.cos(i) * m.cos(aop) - m.cos(raan) * m.sin(aop), 
            m.cos(raan) * m.cos(i) * m.cos(aop) - m.sin(raan) * m.sin(aop), 
            m.sin(i) * m.cos(aop)]
    row2 = [m.sin(raan) * m.sin(i), 
            -m.cos(raan) * m.sin(i), 
            m.cos(i)]
    return np.array([row0, row1, row2])

# Calculate Eccentric Anomaly
def ecc_anomaly(arr, method, tol=1e-8):
    if method == 'Newton':
        # Newton-Raphson method
        Me, e = arr
        if Me < np.pi / 2.0:
            E0 = Me + e / 2.0
        else:
            E0 = Me - e
        
        for n in range(200):  # Max Number of Iterations
            ratio = (E0 - e * m.sin(E0) - Me) / (1 - e * m.cos(E0))
            if abs(ratio) < tol:
                if n == 0:
                    return E0
                else:
                    return E0 - ratio
            E0 = E0 - ratio
        return False  # Did not converge
    
    elif method == 'tae':
        ta, e = arr
        E = 2 * m.atan(m.sqrt((1 - e) / (1 + e)) * m.tan(ta / 2.0))
        return E
    
    else:
        print('Invalid Method')
        return None

# Takes .txt file containing TLE data and returns COEs and epoch date
def tle2coes(tle_filename, mu=pd.Earth['mu'], return_date=False):
    with open(tle_filename, 'r') as file:
        lines = file.readlines()

    # Separate into Three Lines
    line0 = lines[0].strip()  # Name
    line1 = lines[1].strip().split()
    line2 = lines[2].strip().split()

    # Extract Relevant Data from TLE
    # Epoch
    epoch = line1[3]
    year, month, day, hour = calc_epoch(epoch)
    epoch_date = datetime(year, month, day) + timedelta(hours=hour)

    # Classical Orbital Elements
    # Inclination (radians)
    i = float(line2[2]) * d2r
    # Right Ascension of the Ascending Node (radians)
    raan = float(line2[3]) * d2r
    # Eccentricity
    e = float('.' + line2[4])
    # Argument of Perigee (radians)
    aop = float(line2[5]) * d2r
    # Mean Anomaly (radians)
    Me = float(line2[6]) * d2r
    # Mean Motion (revolutions per day)
    mean_motion = float(line2[7])
    # Time Period (seconds)
    T = 1 / mean_motion * 24 * 3600
    # Semi-Major Axis (km)
    a = (mu * (T / (2 * np.pi))**2)**(1 / 3)
    # Eccentric Anomaly
    E = ecc_anomaly([Me, e], 'Newton')
    # True Anomaly
    ta = true_anomaly([E, e])
    # Magnitude of Radius Vector
    r_mag = a * (1 - e * np.cos(E))

    if return_date:
        return [a, e, i, ta, aop, raan], epoch_date
    else:
        return [a, e, i, ta, aop, raan]

def calc_epoch(epoch):
    # Epoch Year
    year = int('20' + epoch[0:2])
    epoch_parts = epoch[2:].split('.')

    # Epoch Day of Year
    day_of_year = int(epoch_parts[0]) - 1

    # Decimal Hour of Day
    decimal_hour = float('0.' + epoch_parts[1]) * 24.0

    # Year Month Day
    date_obj = datetime(year, 1, 1) + timedelta(days=day_of_year)
    month = int(date_obj.month)
    day = int(date_obj.day)

    return year, month, day, decimal_hour

def true_anomaly(arr):
    E, e = arr
    return 2 * m.atan(m.sqrt((1 + e) / (1 - e)) * m.tan(E / 2.0))

def tle2rv(tle_filename):
    return coes2rv(tle2coes(tle_filename))