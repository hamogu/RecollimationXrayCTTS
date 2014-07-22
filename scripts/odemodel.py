'''
We perform the actual fit with the Sherpa fitting tool 
(http://cxc.harvard.edu/sherpa4.6/).
Sherpa is distributed as a complied binary with just a small set of python
modules. It includes numpy, but not scipy. Thus, we cannot run the 
stellarwindshock model in sherpa. Unfortunately, it is difficult to compile
sherpa in a different environment. Thus, this module implements the 
interface required for a sherpa model, but the actual computation of the 
stellar wind shock are done in a separate shell launched in a subprocess.
The PATH in that shell points to a python installation that includes scipy
(the module level variable PATH_TO_PYTHON_WITH_SCIPY contains the PATH for
that environment).
For all computations, this model calls ``connectsherpa.py`` which is a 
driver script, that allows us to run a stellar wind shock model on the 
command line and pass the results back as a pickeled output string.

The driver script also does some trivial computations (e.g. binning the output
to reduce the amount of data that needs to be passed through the pickeled 
stream.

How use this model?
-------------------
See ``model_fit_in_sherpa`` for a usage example.

Essentially, loading this module in Sherpa adds a model ``recol``, that can 
be used as usual:
>>> set_source(1, xsphabs.a1 * recol.r)

To fit the distance between the star and the position of the x-ray shock,
make a dataset with x-value 123456. This is the sign for the model to output 
the distance in AU as y-value. E.g. the following code fits model to an X-ray
position 30+-5 AU from the star:

>>> load_arrays(2,[123456], [30], [5])
>>> set_source(2, recol.r1)


'''


import cPickle
import subprocess

import numpy as np

from sherpa import models
from sherpa.astro import ui

PATH_TO_PYTHON_WITH_SCIPY = '/data/guenther/anaconda/bin'

def model_mass_flux(v, n, d=140.):
    '''calculate the mass flux for a shock model with norm=1

    Shock models (originally for accretion shocks) from Guenther et al, 
    A&A 446, 1111 (2007) updated in Guenther, AN 332, 448 (2011) are
    available from:
    http://hdl.handle.net/10904/10202
    These shock models are normalized to a certain shock area.
    This function calculates the mass flux for a normalized shock in those
    models.

    Parameters
    ----------
    v : float
        pre-shock velocity in km/s
    n : float
        pre-shock density in cm^{-3}
    d : float
        distance to target in pc

    Returns
    -------
    massflux : float
        mass flux of a model with norm=1 in cgs units
    '''
    return 1.67e-24 * n * v*1e5 * 1e20* (d/10.)**2


class Recol(models.ArithmeticModel):
    
    def __init__(self, name='Recol'):
        # convert to cm/s
        self.v = models.Parameter(name, 'v', 750., min=400, max=800, 
                        hard_min=0, units='km/s')
        # maybe use M_sun/yr here and convert
        self.mdot = models.Parameter(name, 'mdot', 1e-9, units='M_sun/yr', min=0, hard_min=0)
        self.p0 = models.Parameter(name, 'p0', 1e-5, units='barye', min=0,  hard_min=0)
        self.p_scale = models.Parameter(name, 'p_scale', 5, units='AU', frozen=True, min=0, hard_min=0)
        self.p_inf = models.Parameter(name, 'p_inf', 1e-7, units='barye', min=0, hard_min=0)
        ui.load_table_model('shock', 'shocks.fits')
        self.shock = ui.get_model_component('shock')
        self.shock.n0 = 1e10

        models.ArithmeticModel.__init__(self, name, 
               (self.v, self.mdot, self.p0, self.p_scale, self.p_inf))


    #@modelCacher1d
    def calc(self, pars, *args):
        mdot = self.mdot.val /  (365.25 * 24. *3600. / 2e33)
        v = self.v.val * 1e5
        pars_in = "{0} {1} {2} {3} {4}".format(v, mdot, self.p0.val, 
                                               self.p_scale.val, self.p_inf.val)
        try:
            out = subprocess.check_output(['python connectsherpa.py '+pars_in],
                             shell=True, 
                             env={'PATH': PATH_TO_PYTHON_WITH_SCIPY})
        except subprocess.CalledProcessError, ValueError:
            print "ValueError"
            return np.zeros_like(args[0])
        oderesult = cPickle.loads(out)
        print oderesult
        # x = 123456 is a special value and means: output z_max
        if np.abs(args[0][0] - 123456) < 0.1:  #allow for rounding errors
            return np.array([oderesult['z_max']])
        else:
              flux = np.zeros_like(args[0])
            for v, m in zip(oderesult['v'], oderesult['massflux']):
                if m > 0: 
                    self.shock.v0 = v
                    self.shock.norm = m/model_mass_flux(v, 1e10)
                    #print shock.v0.val, shock.norm.val
                    flux += self.shock(*args)
            return flux

ui.add_model(Recol)
