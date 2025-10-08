import spiceypy as spice
import numpy as np

# Retrieves IDs, Names, and Time Coverages for all Objects in spk File
def get_objects(filename, display = False):
    objects = spice.spkobj(filename)
    ids, names, tcs_sec, tcs_cal = [], [], [], []
    n = 0
    if display:
        print('\nObjects in %s:' % filename)

    for o in objects:
        ids.append(o)

        # Time Coverage in Seconds since J2000
        tc_sec = spice.spkcov(filename, o)

        # Convert Time Coverage to Calendar Dates
        tc_cal = [spice.timout(f, 'YYYY MON DD HR:MN:SC.### (TDB) :: TDB') for f in tc_sec]

        # Append Time Coverages to Output Lists
        tcs_sec.append(tc_sec)
        tcs_cal.append(tc_cal)

        # Get Object Name from ID
        try:
            names.append(id2body(o))

        except:
            names.append('Unknown')

        if display:
            print('ID: %i\tname: %s\t\ttc: %s ---> %s' % (ids[-1], names[-1], tc_cal[0], tc_cal[1]))
           
    return ids, names, tcs_sec, tcs_cal


# Returns Name of Body from ID
def id2body(id_):
    return spice.bodc2n(id_)

# Creates Time Array for Given Time Coverage
def tc2array(tcs, steps):
    arr = np.zeros((steps, 1))
    arr[:,0] = np.linspace(tcs[0], tcs[1], steps)
    return arr

# Get Ephemeris from a Given Time Array
def get_ephemeris_data(target, times, frame, observer):
    return np.array(spice.spkezr(target, times, frame, 'NONE', observer)[0])

