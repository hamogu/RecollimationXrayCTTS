'''Driver script to calculate a stellar wind shock from the command line

see ``odemodel.py`` for how to use this with the sherpa fitting tool
'''

import os
import sys
import cPickle

import numpy as np

from stellarwindshock import StellarWindShock, ExpPres

# from http://stackoverflow.com/questions/6796492/python-temporarily-redirect-stdout-stderr
class RedirectStdStreams(object):
    def __init__(self, stdout=None, stderr=None):
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr

    
if len(sys.argv) != 6:
    raise ValueError('Calling sequence is > connectsherpa.py v mdot p0 scale p_inf')

v, mdot, p0, scale, p_inf = sys.argv[1:]
v = float(v)
mdot = float(mdot)
p0 = float(p0)
scale = float(scale)
p_inf = float(p_inf)

# hardcode some parameters here to reduce the amount of data that needs to
# pass through pickeling
z = np.arange(0, 10000, .01)
omega_0 = 0.01
vbins = np.arange(250., 1051., 100.)
z_min = 5

pressure = ExpPres(p0, scale, p_inf)
ode = StellarWindShock(pressure, v, mdot)

devnull = open(os.devnull, 'w')
with RedirectStdStreams(stdout=devnull, stderr=devnull):
    # ode.solve will also print convergence problems.
    # need to ignore those here, so only pickeled output gets written
    ode.solve(omega_0, z)

ind = ode.z > z_min
hist, bin_edges = np.histogram(ode.v_shock[ind], bins=vbins*1e5, weights=ode.massflux[ind])
v = (bin_edges[0:-1] + bin_edges[1:])/2.

print cPickle.dumps({'z_max': ode.z.max(), 'v': v/1e5, 'massflux': hist}, 0)
