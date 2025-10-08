# Tools Script

import numpy as np
import matplotlib.pyplot as plt
import planetary_data as pd
import math as m
from datetime import datetime, timedelta

d2r = np.pi / 180.0  # Degrees to Radians
r2d = 180.0 / np.pi  # Radians to Degrees

def norm(v):
    return np.linalg.norm(v)

def plot_n_orbits(rs, labels, cb=pd.Earth):
    ax = plt.axes(projection='3d')
    n = 0
    for r in rs:
        ax.plot(r[:, 0], r[:, 1], r[:, 2], label=labels[n])
        ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], 'o')
        n += 1
    
    r_plot = cb['radius']
    # Create a sphere for Earth
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = r_plot * np.cos(_u) * np.sin(_v)
    _y = r_plot * np.sin(_u) * np.sin(_v)
    _z = r_plot * np.cos(_v)
    ax.plot_surface(_x, _y, _z, color='b', alpha=0.3)
    
    l = r_plot * 2.0
    x, y, z = [0, 0, 0], [0, 0, 0], [0, 0, 0]
    u, v, w = [[l], [0], [0]], [[0], [l], [0]], [[0], [0], [l]]
    ax.quiver(x, y, z, u, v, w, color='r', arrow_length_ratio=0.1)
    
    max_value = np.max(np.abs(rs))
    ax.set_xlim([-max_value, max_value])
    ax.set_ylim([-max_value, max_value])
    ax.set_zlim([-max_value, max_value])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('equal')
    plt.legend()
    plt.show()

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