"""
Figure 7: Gao et al. atmospheric structure

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
    savetag = "fig07"
    fname = "1barGao_final3.atm"
    plot_title = "1 bar CO$_2$/CO/O$_2$, Desiccated"
    seed = 2
    xlim = [6e-8, 3]
    ylim = [1e5, 1e-2]
    tlim = [180, 300]

    # Read in atm files
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname)
    data = np.genfromtxt(atmpath, skip_header=2)

    # Parse atm data
    P, T, alt = data[:,0], data[:,1], data[:,2]
    gas_profiles = data[:,3:]
    molec_names = ["O$_2$", "O$_3$", "CO", "CO$_2$", "H$_2$O", "H$_2$O$_2$", "OH", "H$_2$", "HO$_2^-$"]

    # Feed to plotting function
    #atm.add_atm_plot(P, T, gas_profiles, molec_names, legend=True, title=plot_title)
    atm.plot_single_atm(P, T, gas_profiles, molec_names, title=plot_title,
                        savetag=savetag, seed=seed, xlim=xlim, ylim=ylim, tlim=tlim)

    return
#########################################

if __name__ == "__main__":

    make_fig()
