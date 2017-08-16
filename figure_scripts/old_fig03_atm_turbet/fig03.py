"""
Figure 3

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
    savetag = "fig03"
    #fname = "figure3.txt"
    fname = "turbet_fixed_h2o.atm"
    plot_title = "1 bar N$_2$, 376 ppm CO$_2$, Moist"
    seed = 2
    xlim = [1e-9, 3]
    ylim = [1e5, 1e-2]
    tlim = [155, 325]

    # Read in atm files
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname)
    data = np.genfromtxt(atmpath, skip_header=2)

    # Parse atm data
    P, T, alt = data[:,0], data[:,1], data[:,2]
    gas_profiles = data[:,3:]
    # set N2 to fill
    gas_profiles[:,-1] = 1e0 - gas_profiles[:,0] - gas_profiles[:,1]
    molec_names = ["CO$_2$", "H$_2$O", "N$_2$"]

    # Feed to plotting function
    #atm.add_atm_plot(P, T, gas_profiles, molec_names, legend=True, title=plot_title)
    atm.plot_single_atm(P, T, gas_profiles, molec_names, title=plot_title,
                        savetag=savetag, seed=seed, xlim=xlim, ylim=ylim, tlim=tlim)

    return
#########################################

if __name__ == "__main__":

    make_fig()
