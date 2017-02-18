"""
Figure 10

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    import numpy as np
    import os
    import sys
    sys.path.append("../../")
    from utils import fig_params, atm

    # Params specific to this plot
    savetag = "fig10"
    fname = "fig10_hcaer4.out"

    # Read in atm files
    path = os.path.join(os.path.dirname(__file__),"model_outputs/", fname)
    data1 = np.genfromtxt(path, skip_header=4)

    # Parse atm data. cm-->km, , cm-->um
    alt, aero, rpar = data1[:,0]*1e-5, data1[:,1], data1[:,2]*1e4

    # Plot
    atm.plot_aerosol(alt, aero, rpar, savetag=savetag)

    return
#########################################

if __name__ == "__main__":

    make_fig()
