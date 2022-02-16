import numpy as np
import numpy.ctypeslib as npct
import ctypes

array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')

libpiston = npct.load_library('libpiston', '.')

libpiston.simulate_piston.restype = None
libpiston.simulate_piston.argtypes = [array_1d_double, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]

def simulate(timeSteps, deltat, Npart, mpart, T, X0, gravityg):
    positions = np.zeros(timeSteps)
    times = np.arange(timeSteps)*deltat
    libpiston.simulate_piston(positions, timeSteps, deltat, Npart, mpart, T, X0, gravityg)

    return times, positions
