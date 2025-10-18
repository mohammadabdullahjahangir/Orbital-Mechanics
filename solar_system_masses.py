# Masses of Bodies in Solar System

import sys 
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
import tools as t
import spiceypy as spice
import spice_tools as st

# Gravitational Constant
G_meters = 6.67430e-11  # km^3/kg/s^2, gravitational constant
G = G_meters *10**-9

# SPICE Body IDs for Planets (excluding Mercury and Venus)
PLANETS = range(1, 10)

# Colours for Plot
colours = ['C1', 'C2', 'g', 'r', 'C3', 'b', 'c', 'm', 'y', 'C4']

if __name__ == '__main__':

    # Load Kernel
    spice.furnsh('spice_data/gm_de431.tpc')

    # Initialise Figure
    plt.figure(figsize = (10, 5))

    # For each SPICE Body ID
    for n in PLANETS:

        # Get Mass of Planet and Barycenter of System
        mass_planet = spice.bodvcd(n * 100 + 99, 'GM', 1)[1][0]
        mass_barycenter = spice.bodvcd(n, 'GM', 1)[1][0]

        # Calculate Percentage Ratio of Planet to Total System Mass
        ratio = mass_planet / mass_barycenter * 100

        # Print Values to Terminal
        print(spice.bodc2n(n))
        print(ratio)
        print('\n')

        # Add Values to Plot
        plt.semilogx([mass_planet / G], [ratio], colours[n - 1] + 'o')

        # Even Values go Below Point, Add Above
        if n == 6:
            xytext = (0, -15)
        elif n % 2 == 1:
            xytext = (0, -30)
        else:
            xytext = (0, 10)

        # Add Annotation to Plot
        plt.annotate(spice.bodc2s(n * 100 + 99), [mass_planet / G, ratio], textcoords = 'offset points', xytext = xytext, ha = 'center', color = colours[n - 1], fontsize = 'large', fontweight = 'bold')


    # List for all Masses in Solar System
    mass_ss = np.zeros(10)

    # For each Planet ID (including Sun)
    for n in range(1, 11):

        # Add Mass to List
        mass_ss[n - 1] = spice.bodvcd(n, 'GM', 1)[1][0]

    # Ratio of Mass to Sun
    ratio = mass_ss[-1] / sum(mass_ss) * 100

    # Add Point to Plot
    plt.semilogx([mass_ss[-1] / G], [ratio], colours[-2] + 'o')

    # Add Annotation
    plt.annotate('SUN', [mass_ss[-1] / G, ratio], textcoords = 'offset points', xytext = xytext, ha = 'center', color = colours[-2 ], fontsize = 'large', fontweight = 'bold')

    # Parameters
    plt.ylim([87.5, 102])
    plt.title(' Planet System Mass Ratios')
    plt.ylabel('Mass Ratio', fontsize = 12)
    plt.xlabel('Planet Mass (kg)', fontsize = 12)
    plt.grid()
    # plt.savefig('mass_ratios.png', dpi = 400)
    plt.show()