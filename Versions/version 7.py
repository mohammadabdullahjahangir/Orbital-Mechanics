# Using Two Line Element Sets (TLEs) to Track Satellites

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP


# Time Parameters
tspan = 24.0 * 3600.0 
dt = 10.0

if __name__ == '__main__':

    # ISS
    op0 = OP(t.tle2coes('iss.txt'), None, None, tspan, dt, coes = True, deg = False)
    
    # Hubble
    op1 = OP(t.tle2coes('hubble.txt'), None, None, tspan, dt, coes = True, deg = False)

    # Praetorian
    op2 = OP(t.tle2coes('praetorian sda_613.txt'), None, None, tspan, dt, coes = True, deg = False)

    # Starlink
    op3 = OP(t.tle2coes('starlink_35174.txt'), None, None, tspan, dt, coes = True, deg = False)

    t.plot_n_orbits([op0.rs, op1.rs, op2.rs, op3.rs], ['ISS', 'Hubble', 'Praetorian', 'Starlink'])