# Keplerian Orbit Elements

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP

style.use('dark_background')

# Time Parameters
tspan = 24.0 * 3600.0 
dt = 10.0

# Central Body
cb = pd.Earth

if __name__ == '__main__':
    # a, e, i, ta, aop, raan = coes
    # ISS
    c0 = [cb['radius'] + 417.5, 0.0000898, 51.6318, 0.0, 206.6768, 115.2314]

    # Geostationary
    c1 = [cb['radius'] + 35786.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Random
    c2 = [cb['radius'] + 20000.0, 0.5, 30.0, 0.0, 45.0, 60.0]

    op0 = OP(c0, None, None, tspan, dt, coes = True)
    op1 = OP(c1, None, None, tspan, dt, coes = True)
    op2 = OP(c2, None, None, tspan, dt, coes = True)

    op0.propagate_orbit()
    op1.propagate_orbit()
    op2.propagate_orbit()

    t.plot_n_orbits([op0.rs, op1.rs, op2.rs], ['ISS', 'Geostationary', 'Random'])
