"""
This file contains global figure parameters that are used for all plots in
this project.

Author: Jacob Lustig-Yaeger
"""

import matplotlib as mpl

mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rc('axes', linewidth=1.5)
mpl.rc('lines', linewidth=2)
mpl.rcParams['font.size'] = 30.0

COLOR1 = "black"
LEGEND_FONTSIZE = 16

def colors_from_molecs(molecules):
    return

def set_default_figsize():
    mpl.rcParams['figure.figsize'] = (10,8)

def set_atm_figsize():
    mpl.rcParams['figure.figsize'] = (7,8)

def set_spectrum_figsize():
    mpl.rcParams['figure.figsize'] = (16,8)

def tex_molecule(formula):
    """
    Create LaTeX compatible string for molecular formula from simple string
    """
    new = r''
    for char in formula:
        if char.isdigit():
            tmp = '$_{'+char.upper()+'}$'
        else:
            tmp = char
        new = new+tmp
    return new
