"""
This file contains global figure parameters that are used for all plots in
this project.

Author: Jacob Lustig-Yaeger
"""

import matplotlib as mpl

# Typical plot parameters that make for pretty plots
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 30.0

## for Palatino and other serif fonts use:
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rc('text', usetex=True)

color1 = "black"
