# This are the commands that I used to get the data for DG Tau
# They assume that the CIAO system is started already and should be executed
# in a shell (not in Python)

# > download_chandra_obsid 4487,6409,7246,7247

# > specextract "*/*/*evt2*[sky=circle(04:27:04.70,26:06:16.3,0.05')]" DGTau bkgfile="*/*/*evt2*[sky=circle(04:27:07.6,26:09:28,1.4')]" combine=yes

# These are the commands I used to fitting in Sherpa

import odemodel
load_data('/data/guenther/obs/Chandra/DG_Tau/DGTau_combined_src.pi')
group_counts(20)
subtract()
ignore(None, 0.3)
ignore(7., None)
set_source(xsphabs.a1 * xsvapec.v1+ xsphabs.a2 * xsvapec.v2)
# Start close to Guedel et al values, so that jet is comp1 and corona is comp2
a1.nh = 0.1
v1.kT = 0.2
a2.nH = 1
v2.kT = 2
fit(1)
freeze(a2.nh,v2.kT, v2.norm)

load_arrays(2,[123456], [30], [5])
set_source(2, recol.r1)
r1.p_inf=r1.p0 * 0.01
r1.p0.min=1e-7
r1.p0.max=1e-4
r1.mdot.max=1e-6
r1.mdot.min=1e-10
r1.mdot = 5e-10
r1.p_scale=5
r1.v.max=1000
set_source(1, xsphabs.a1 * recol.r1 + xsphabs.a2 * xsvapec.v2)
# set as close to best fit values as I can
a1.nH     =     0.466986    
r1.v      =     837.509     
r1.mdot   =     5.17862e-10 
r1.p0     =     1.43298e-05 


plot_fit()
plot_model_component(a1*r1, overplot=True)
plot_model_component(a2*v2, overplot=True)
log_scale()
limits(Y_AXIS, 0.0001, 0.02)
limits(X_AXIS, .55, 6.4)
set_plot_title("")
printopts = {'bottommargin':0.1, 'topmargin':.02, 'leftmargin':.1, 'rightmargin': .02, 'dpi': 300, 'clobber': True, 'fittopage': True}

printopts = {'scaleheight':500, 'dpi': 300, 'clobber': True, 'fittopage': True}
# because my computer has an ATI Graphics card that turn the plot window ?!?
set_preference("export.printmethod", "xoffscreen")
set_plot( { 'bottommargin' :0.15, 'topmargin': 0.02, 'leftmargin': 0.2, 'rightmargin': 0.02 } )
set_axis("ax1","ticklabel.size=30 label.size=30 offset.perpendicular=50.00")
set_axis("ay1","ticklabel.size=30 label.size=30 offset.perpendicular=100.00")

print_window('../figures/DGTaufit/DGTaufit.png', {'dpi': 300, 'clobber': True})
print_window('../figures/DGTaufit/DGTaufit.eps', printopts)


"""
Below, I paste the ourput of the best fit model and the conf()
command in sherpa after finting the best fit.
These values appear in the table and in the text in the article.

{'z_max': 30.350000000000001, 'massflux': array([  5.50812136e+14,   6.20901988e+15,   2.87946247e+15,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00]), 'v': array([  300.,   400.,   500.,   600.,   700.,   800.,   900.,  1000.])}
Datasets              = 1, 2
Method                = levmar
Statistic             = chi2gehrels
Initial fit statistic = 290.308
Final fit statistic   = 52.7545 at function evaluation 176
Data points           = 50
Degrees of freedom    = 46
Probability [Q-value] = 0.229291
Reduced statistic     = 1.14684
Change in statistic   = 237.554
   a1.nH          0.466986    
   r1.v           837.509     
   r1.mdot        5.17862e-10 
   r1.p0          1.43298e-05 
Model: 1
((xsphabs.a1 * recol.r1) + (xsphabs.a2 * xsvapec.v2))
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   a1.nH        thawed     0.466986            0       100000 10^22 atoms / cm^2
   r1.v         thawed      837.509          400         1000       km/s
   r1.mdot      thawed  5.17862e-10        1e-10        1e-06   M_sun/yr
   r1.p0        thawed  1.43298e-05        1e-07       0.0001      barye
   r1.p_scale   frozen            5            0  3.40282e+38         AU
   r1.p_inf     linked  1.43298e-07     expr: (r1.p0 * 0.01)      barye
   a2.nH        frozen      2.61658            0       100000 10^22 atoms / cm^2
   v2.kT        frozen      2.31504       0.0808       68.447        keV
   v2.norm      frozen  0.000203099            0        1e+24 
r1.p0 +: WARNING: The confidence level lies within (1.543357e-05, 1.543270e-05)
r1.p0 upper bound:	1.10334e-06
Datasets              = 1, 2
Confidence Method     = confidence
Iterative Fit Method  = None
Fitting Method        = levmar
Statistic             = chi2gehrels
confidence 1-sigma (68.2689%) bounds:
   Param            Best-Fit  Lower Bound  Upper Bound
   -----            --------  -----------  -----------
   a1.nH            0.466986   -0.0262755    0.0279926
   r1.v              837.509     -7.15593      22.1393
   r1.mdot       5.17862e-10 -1.60783e-10  3.23427e-11
   r1.p0         1.43298e-05 -6.71709e-07  1.10334e-06

"""
