import os

# Get the base directory of the project
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Define paths to SPICE kernels
leap_seconds_kernel = os.path.join(BASE_DIR, 'spice_data', 'latest_leapseconds.tls')
de432s = os.path.join(BASE_DIR, 'spice_data', 'de432s.bsp')