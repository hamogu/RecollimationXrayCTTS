'''Recollimation boundary layers as X-ray sources in young stellar jets

This module contains Python code to calculate the properties of recollimation
shocks of a stellar wind that are caused by an external pressure (e.g. a
disk wind) as described in the manuscript itself.
This module provides classes for different pressure profiles to be used with
either a hot or a cold stellar wind.
Individual functions are documented, for an example see the IPython notebook
called `jet_recol_shocks.ipynb`.
'''
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d

try:
    from matplotlib import pyplot as plt
    _has_mpl = True
except ImportError:
    _has_mpl = False


class NoMatplotlibException(Exception):
    pass


class ConstantPres(object):
    '''A constant pressure'''
    def __init__(self, p_const):
        '''
        Parameters
        ----------        
        p_inf : float
            presure
        '''
        self.p_const = p_const
       
    def __call__(self, z):
       '''
       Parameters
       ----------
       z : array-like
           same dimension as z_1
       '''
       return self.p_const + np.zeros_like(z)


class PowerPres(object):
    '''A powerlaw pressure profile
   
    p(z) = p_inf + p_1 * ((z_1 + z)/z_1)**(-eta)
    '''
    def __init__(self, p_inf, p_1, eta, z_1):
        '''
        Parameters
        ----------        
        p_inf, p_1, eta, z_1 : float
        '''
        self.p_inf = p_inf
        self.p_1 = p_1
        self.eta = eta
        self.z_1 = z_1
       
    def __call__(self, z):
       '''
       Parameters
       ----------
       z : array-like
           same dimension as z_1
       '''
       return self.p_inf + self.p_1 * ((self.z_1 + z)/self.z_1)**(-self.eta)


class ExpPres(object):
    '''An exponential pressure profile

    p_inf + p_0 * np.exp(-z/h)
    '''
    def __init__(self, p_0, h, p_inf = 0.):
        '''
        p(z) = p_inf + self.p_0 * np.exp(-z/h)

        Parameters
        ----------        
        p_0, h, p_inf : float
        '''
        self.p_0 = p_0
        self.h = h
        self.p_inf = p_inf
       
    def __call__(self, z):
        '''
        Parameters
        ----------
        z : array-like
            same dimension as z_1
        '''
        return self.p_inf + self.p_0 * np.exp(-z/self.h)



'''x is pre-shock velocity in km/s, result is fraction of kinetic energy 
radiated in X-rays. The values for the y-axis (the fraction of the total 
luminosity emitted in X-rays) are justified in the appendix of the notebook.
Here they are just copied verbatim, to ensure that this code runs even if the 
shock models of Guenther et al. (2007) are not installed.
'''
xrayfrac = interp1d(np.arange(0.,1001.,100.),
                    np.array([0, 0, 0, 0.02, 0.10, 0.24, 0.40, 0.51, 0.57, 0.61, 0.65]))


class ODEnotsolvedException(Exception):
    '''ODE needs to be solved before properties of the solution are available'''
    pass


class StellarWindShock(object):
    '''This class describes a shock front in a confinement by presure.

    An object of this class describes one specific combination of stellar wind
    parameters and external presure profile. This class provides methods to 
    solve the ordinary differential equation that governs the shape of the 
    shock front, which forms due to the external presure from the expanding 
    disk wind.

    ..note::
      This class accepts all lengths scales input in units of AU and all other
      parameters (velocity, mass flux, presure) in cgs units.
      The possibility for user error would be greatly decreased by using a
      module that tracks the units (e.g. ``astropy.units``) but the scipy ODE
      solver does not work well with any such scheme. Thus, BEWARE of the UNITS!
    '''
    AU = 1.5e13
    k_B = 1.38e-16  # Boltzman constant
    m_H = 1.67e-24  # mass of hydrogen
    mu = 0.7    # mean particle mass in highly ionized plamsa
    def __init__(self, p_ext, v_inf, m_dot):
        '''
        Parameters
        ----------
        p_ext: function or callable object
            ``p_ext`` must accept a height (in AU) and return a pressure 
            (in cgs units)
        v_inf: float
            wind speed (in cgs units)
            It is assumed that the stellar wind reaches its final wind speed 
            before any interaction with the disk wind occurs.
        m_dot: float
            Spherically integrated stellar mass loss rate (in cgs units)
        '''
        self.p_ext = p_ext
        self.v_inf = v_inf
        self.m_dot = m_dot
   

    @property
    def m_dot_sun_yr(self):
        return self.m_dot * 3600.*24.*365.25/2e33


    def R0(self, z):
        '''Eqn 14'''
        const = self.m_dot * self.v_inf / (4.*np.pi* self.p_ext(z))
        # const from cgs in AU         
        return np.sqrt(const)/self.AU
      

    def __call__(self, omega, z):
        r = np.sqrt(omega**2+z**2)
        R0 = self.R0(z)
        r = np.clip(r, 0, R0)
        # eqn. 13
        return np.tan(np.arctan2(omega, z) - np.arcsin(r/R0))
   

    def start(self):
        '''gives equilibrium solution at z=0
   
        Returns
        -------
        d : float
            omega (in AU) for equilibrium solution at z=0
        '''
        return self.R0(0.)

   
    def solve(self, omega0, z):
        '''solve ODE starting at omega=omega0 for z = z[0]

        Parameters
        ----------
        omega0: float
            distance from star (in AU) along the radial direction of the
            cylindrical coordinate system where to start the integration
        z: np.array
            grid in z direction (in AU)
        '''
        if omega0 > self.start():
            s = 'omega0 {0:5.2f} is outside the physically possible range {1:5.2f}'
            raise ValueError(s.format(omega0, self.start()))
               
        self.z = z
        self.omega0 = omega0
        res, info = odeint(self, omega0, z, full_output = True)
        self.res = res.flatten()
        self.odeinfo = info

        # This is to cut solution after shock merged z-axis
        if np.any(self.res < 0.001):
            ngood = np.min(np.where(self.res<0.001)[0])
            self.res = self.res[0:ngood]
            self.z = self.z[0:ngood]
   

    def stellar_wind_density(self, omega, z):
        '''calculate the density of an undisturbed stellar

        see eqn 6 in paper
       
        Parameters
        ----------
        omega: float or np.array
            distance from star (in AU) along the radial direction of the
            cylindrical coordinate system
        z: float or np.array
            distance z direction (in AU)
       
        Returns
        -------
        rho: float or np.array
            density (in cgs units) of the stellar wind
        '''
        if len(omega) != len(z):
            raise ValueError('omega and z must have the same number of elements')
        r = np.sqrt(omega**2 + z**2) * self.AU
        return self.m_dot / (4.*np.pi * r**2. * self.v_inf)
       
  
    @property
    def v_shock(self):
        '''pre-shock wind velocity component perpendicular to shock front'''
        try:
            v_shock = self.v_inf * np.sin(np.arctan2(self.res, self.z)
                                          - np.arctan(self(self.res, self.z)))
        except AttributeError:
            raise ODEnotsolvedException("ODE needs to be solved first")
        return v_shock
   

    @property
    def T_postshock(self):
        '''Post-shock temperature in K (eqn 15)'''
        return 3./16.*self.mu*self.m_H/self.k_B*self.v_shock**2


    @property
    def massflux(self):
        '''numerical approximation to mass flux for each step in the solution'''
        try:
            theta = np.arctan2(self.res, self.z)
        except AttributeError:
            raise ODEnotsolvedException("ODE needs to be solved first")
        # This is calculated by bin, while really it's a continuous quantity.
        # Using np.diff() results in a vector one element short.
        delta_theta = np.zeros_like(theta)
        # theta is reverse ordered, so np.diff is negative
        delta_theta[:-1] = np.abs(np.diff(theta))
        return self.m_dot * np.sin(theta) * delta_theta
   

    @property
    def L_X(self):
        '''X-ray luminosity in erg/s for each step in z'''
        return 0.5 * self.massflux * self.v_shock**2*xrayfrac(self.v_shock/1e5)


    @staticmethod
    def setup_plot(figsize=(8.5,2)):
        '''set up 4 panels to plot some parameters of an ODE solution

        Returns
        -------
        fig: figure instance
        axes: list of for axis instances
        '''
        if _has_mpl:
            fig = plt.figure(figsize=figsize)
            ax1 = fig.add_axes([.06,.24,.21,.7])
            ax2 = fig.add_axes([.27,.24,.21,.7], sharey = ax1)
            ax3 = fig.add_axes([.48,.24,.21,.7], sharey = ax1)
            ax4 = fig.add_axes([.77,.24,.21,.7])
            plt.setp([a.get_yticklabels() for a in fig.axes[1:-1]], visible=False)

            #ax1.set_title('External pressure')
            ax1.set_ylabel('Height (AU)')
            ax1.set_xlabel(r'$P_{\mathrm{ext}}$ (Ba)')
            ax1.set_xscale('log')

            #ax2.set_title('Shock location')
            ax2.set_xlabel('shock position $\omega$ (AU)')

            #ax3.set_title('Shock velocity')
            ax3.set_xlabel('$v_0$ (km s$^{-1}$)')

            #ax4.set_title('Post-Shock temperature')
            ax4.set_xlabel(r'$T_{\mathrm{post-shock}}$ (MK)')
            ax4.set_ylabel('mass fraction per bin')

            return fig, [ax1, ax2, ax3, ax4]
        else:
            raise NoMatplotlibException('Plotting requires Matplotlib')
   

    def plot_params(self, axes, bins=np.arange(0,10.,1.), **kwargs):
        '''plot P_ext, position, v_shock and a T_posthock histrogram in 4 panels

        Parameters
        ----------
        axes : list of four axis instances
        bins : see plt.hist()
               To compare different histograms correctly, make sure to set the 
               same bin boundaries for all of them.
       
        Additional parameters will be passed to ``plt.plot`` and ``plt.hist``.
        Use this for simple customization of plots (e.g. colors).
        '''
        if _has_mpl:
            axes[0].plot(self.p_ext(self.z), self.z, **kwargs)
            axes[1].plot(self.res, self.z, **kwargs)
            axes[2].plot(self.v_shock/1e5, self.z, **kwargs)
            axes[3].hist(self.T_postshock/1e6, 
                         weights=self.massflux/self.massflux.sum(),
                         histtype='step', bins=bins, **kwargs)
        else:
            raise NoMatplotlibException('Plotting requires Matplotlib')


    def plot(self):
        '''plot P_ext, position, v_shock and a T_postshock histogram in 4 panels
        '''
        if _has_mpl:
            fig, axes = self.setup_plot()
            self.plot_params(axes)
            axes[1].plot(np.sqrt(self.R0(self.z)**2-self.z**2), self.z,
                         'k', label='maximum radius')
            axes[2].set_xlim([0, None])
            axes[3].set_ylim([0, None])
        else:
            raise NoMatplotlibException('Plotting requires Matplotlib')


class HotStellarWindShock(StellarWindShock):
    '''This class extends class:`StellarWindShock` for hot winds.
   
    This class is derived from class:`StellarWindShock`. It is implemented as a
    separate class rather than just adding a parameter for the initial wind 
    temperature to class:`StellarWindShock` because there are situations when
    no numerical solution is possible within the assumptions that went into 
    the implemented formulas.
   
    See discussion in the paper for details.
    '''
    def __init__(self, p_ext, v_inf, m_dot, T_0):
        '''
        Parameters
        ----------
        p_ext: function or callable object
            ``p_ext`` must accept a height (in AU) and return a pressure 
            (in cgs units)
        v_inf: float
            wind speed (in cgs units)
            It is assumed that the stellar wind reaches its final wind speed 
            before any interaction with the disk wind occurs.
        m_dot: float
            Spherically integrated stellar mass loss rate (in cgs units)
        T_0 : float
            initial wind temperature (in K)
        '''
        super(HotStellarWindShock, self).__init__(p_ext, v_inf, m_dot)
        self.T_0 = T_0
  

    def __call__(self, omega, z):
        r = np.sqrt(omega**2+z**2)
        R0 = self.R0(z)
        r = np.clip(r, 0, R0)
        last_term = np.sqrt(r**2/R0**2-self.k_B*self.T_0/self.mu/self.m_H/self.v_inf**2)
        return np.tan(np.arctan2(omega, z) - np.arcsin(last_term))
   

    def start(self):
        '''gives equilibrium solution at z=0
   
        Returns
        -------
        d : float
            omega (in AU) for equilibrium solution at z=0
        '''
        term2 = self.v_inf**2. + self.k_B * self.T_0 / (self.mu * self.m_H)
        const = self.m_dot / (self.v_inf * 4.*np.pi* self.p_ext(0)) * term2
        # const from cgs in AU         
        return np.sqrt(const)/1.5e13


    @property
    def T_postshock(self):
        '''Post-shock temperature in K'''
        return self.T_0 + 3./16.*0.7*1.67e-24/1.38e-16*self.v_shock**2

    def plot(self):
        '''plot P_ext, position, v_shock and a T_posthock histrogram in 4 panels
        '''
        
        if _has_mpl:
            fig, axes = self.setup_plot()
            self.plot_params(axes)
            z = self.z[self.z <= self.R0(self.z)]
            axes[1].plot(np.sqrt(self.R0(z)**2-z**2), z, 'k', label='maximum radius')
            axes[2].set_xlim([0, None])
            axes[3].set_ylim([0, None])
        else:
            raise NoMatplotlibException('Plotting requires Matplotlib')
