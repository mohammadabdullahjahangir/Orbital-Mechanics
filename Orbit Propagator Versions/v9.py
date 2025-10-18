# Converting State Vector to Keplerian Elements (Helper for Plotting COEs)

import sys 
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

# Time Parameters
tspan = 24 * 3600 * 20
dt = 100.0

cb = pd.Earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2'] = True
    op = OP(t.tle2coes('txt/iss.txt'), None, None, tspan, dt, coes = True, perts = perts)
    op.calculate_coes()