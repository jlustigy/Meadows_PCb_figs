"""
Figure 4

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
    savetag = "fig04"
    fname1 = "profile_o2lb_10bar_h2o.pt_filtered.atm"
    fname2 = "profile_o2lb_10bar_dry.pt_filtered.atm"
    title1 = "10 bar O$_2$, Ocean"
    title2 = "10 bar O$_2$, Dessicated"
    seed = 7
    xlim = [1e-9, 1.1]
    ylim = [1e6, 1e-3]
    tlim = [155, 325]

    # Read in atm files
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname1)
    data1 = np.genfromtxt(atmpath, skip_header=1)
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname2)
    data2 = np.genfromtxt(atmpath, skip_header=1)

    # Parse atm data
    P1, T1 = data1[:,0], data1[:,1]
    gas_profiles1 = data1[:,2:]
    molec_names1 = ["CO$_2$", "O$_2$", "O$_3$", "CO", "HNO$_3$", "NO$_2$", "SO$_2$", "N$_2$"]

    P2, T2 = data2[:,0], data2[:,1]
    gas_profiles2 = data2[:,2:]
    molec_names2 = ["CO$_2$", "O$_2$", "O$_3$", "CO", "HNO$_3$", "NO$_2$", "SO$_2$", "N$_2$"]

    P = (P1, P2)
    T = (T1, T2)
    gas_profiles = (gas_profiles1, gas_profiles2)
    molec_names = (molec_names1, molec_names2)
    title = (title1, title2)

    # Feed to plotting function
    #atm.add_atm_plot(P, T, gas_profiles, molec_names, legend=True, title=plot_title)
    atm.plot_double_atm(P, T, gas_profiles, molec_names, title=title,
                        savetag=savetag, xlim=xlim,
                        ylim=ylim, tlim=tlim, seed=seed)

    return
#########################################

if __name__ == "__main__":

    make_fig()
