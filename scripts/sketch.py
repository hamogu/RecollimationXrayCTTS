'''Make Fig 1 in the paper

This figure could have been drawn in a graphics program, but we found it 
more convenient to use a short python script.
''' 

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches



fig = plt.figure(figsize=(5,6))
ax = fig.add_subplot(111)
ax.arrow(0.05,0.05,0,.8, lw=4, head_width=0.05, fc='k', length_includes_head=True)
ax.arrow(0.05,0.05,.8,0, lw=4, head_width=0.05, fc='k', length_includes_head=True)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
for i in ['top', 'bottom', 'right', 'left']:
    ax.spines[i].set_visible(False)
ax.text(0.05, .9, '$z$', fontsize=30, horizontalalignment='center', verticalalignment='center')
ax.text(0.9, .05, r'$\omega$', fontsize=30, horizontalalignment='center', verticalalignment='center')

def shock_surface(x):
    '''calculate a function that makes it easy to define the variables of the geometry.'''
    return 4*x**2

x = np.arange(0,0.45,0.01)
ax.plot(0.05+0.15+x, 0.05+shock_surface(x), lw=3, c='k')
ax.text(0.65, .9, r'$R_{\rm{shock}}(z)$', fontsize=30, horizontalalignment='center', verticalalignment='center')

ax.arrow(0.05,0.05,.15+.3, shock_surface(0.3), lw=2, head_width=0.02, fc='k', length_includes_head=True)

ax.plot([.5, 0.5], [.1,.6], ':k', lw=2)
x = np.array([-.1,.1])
ax.plot([.375,.665],[.1,.8], ':k', lw=2)

ax.text(0.1, .15, r'$\theta$', fontsize=30, horizontalalignment='center', verticalalignment='center')
path = Path([(.05,.25),(.12,.25),(.17,.15)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
patch = patches.PathPatch(path, facecolor='none', lw=2)
ax.add_patch(patch)

ax.text(0.47, .25, r'$\alpha$', fontsize=30, horizontalalignment='center', verticalalignment='center')
path = Path([(.5,.2),(.46,.18),(.42,.2)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
patch = patches.PathPatch(path, facecolor='none', lw=2)
ax.add_patch(patch)

ax.text(0.36, .22, r'$\psi$', fontsize=30, horizontalalignment='center', verticalalignment='center')
path = Path([(.4,.17),(.3,.15),(.25,.21)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
patch = patches.PathPatch(path, facecolor='none', lw=2)
ax.add_patch(patch)

x = np.arange(0,0.75,0.01)
ax.plot(0.05+0.15+x, 0.05+1.2*x**2, lw=3, c='k')
ax.text(0.87, .75, r'$R_{\rm{contdisc}}(z)$', fontsize=20, horizontalalignment='center', verticalalignment='center')

ax.text(.25,.7, 'stellar\nwind\n(pre-shock)', fontsize=20, ha='center')
ax.text(.7,.5, 'stellar\nwind\n(post-\nshock)', fontsize=16, ha='center')
ax.text(.75,.15, 'disk\nwind', fontsize=20, ha='center')

ax.set_xlim([0,1])
ax.set_ylim([0,1])
fig.subplots_adjust(left=0.,right=1.,bottom=0,top=1)
fig.savefig('../figures/sketch/sketch.png')
fig.savefig('../figures/sketch/sketch.eps')

