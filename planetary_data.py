import numpy as np

# Converting Days to Seconds
day2sec = 24 * 3600

G_meters = 6.67408e-11  # km^3/kg/s^2, gravitational constant
G = G_meters *10**-9

Sun = {
    'name': 'Sun',
    'radius': 697500.0,
    'mass': 1.981e30,
    'mu': 1.981e30*G,
    'G1': 10.0**8, #kg*km**3/s**2m/**2,
    'spice_file': 'spice_data/de432s.bsp',
    'deorbit_altitude': 2*695510.0,
    # 'J2': 1.081874e-3
}


atm = np.array([[63.096, 2.059e-4], [251.189, 5.909e-11], [1000.0, 3.561e-15]])
Earth = {
    'name': 'Earth',
    'radius': 6378.0,
    'mass': 5.972e24,
    'mu': 5.972e24*G,
    'J2': 1.082635854e-3,
    'zs': atm[:,0],
    'rhos': atm[:,1]*10**8,
    'atm_rot_vector': np.array([0, 0, 7.2921150e-5]),  # rad/s
    'deorbit_altitude': 10.0,
    'spice_file': 'spice_data/de432s.bsp'
}

Moon = {
    'name': 'Moon',
    'radius': 1737.1,
    'mass': 7.34767309e22,
    'mu': 7.34767309e22*G,
    'orbit_T': (29 * day2sec) + (12 * 3600) + (44 * 60) + 2.9,
    'spice_file': 'spice_data/de432s.bsp'
}

Moon['orbit_w'] = 2 * np.pi / Moon['orbit_T']
