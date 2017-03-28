"""
Figure 8

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
    savetag = "fig08_new"
    #fname1 = "profile_modern_earth_preindustrial.pt_filtered.atm"
    #fname2 = "profile_earth_prox.pt_filtered.atm"
    fname1 = "profile_Earth_modern_preindustrial_R1.pt_filtered.atm"
    fname2 = "profile_Earth_proxb_.pt_filtered_R1.atm"
    title1 = "Pre-Industrial Earth"
    title2 = "Earth-like Proxima Cen b"
    seed = 25#16#11#7
    xlim = [1e-8, 2]
    ylim = [1e5, 4e-1 ]
    tlim = [145, 305]
    legloc = [4e3, 1e1, 1e2, 2e2, 5, 1e1, 1, None, None, None, None, 4e2, 3]

    # Read in atm files
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname1)
    data1 = np.genfromtxt(atmpath, skip_header=1)
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname2)
    data2 = np.genfromtxt(atmpath, skip_header=1)

    # Parse atm data
    P1, T1 = data1[:,0], data1[:,1]
    gas_profiles1 = data1[:,2:]
    molec_names1 = ["H$_2$O", "CH$_4$", "CO$_2$", "O$_2$", "O$_3$", "CO", "HNO$_3$", "NO$_2$", "SO$_2$", "N$_2$O", "N$_2$"]
    molec_names1 = ["H$_2$O", "CH$_4$", "CO$_2$", "O$_2$", "O$_3$", "CO", "H$_2$CO", "HNO$_3$", "NO$_2$", "SO$_2$", "N$_2$O", "N$_2$"]

    P2, T2 = data2[:,0], data2[:,1]
    gas_profiles2 = data2[:,2:]
    molec_names2 = ["H$_2$O", "CH$_4$", "CO$_2$", "O$_2$", "O$_3$", "CO", "HNO$_3$", "NO$_2$", "SO$_2$", "N$_2$O", "N$_2$"]
    molec_names2 = ["H$_2$O", "CH$_4$", "CO$_2$", "O$_2$", "O$_3$", "CO", "H$_2$CO", "HNO$_3$", "NO$_2$", "SO$_2$", "N$_2$O", "N$_2$"]

    P = (P1, P2)
    T = (T1, T2)
    gas_profiles = (gas_profiles1, gas_profiles2)
    molec_names = (molec_names1, molec_names2)
    title = (title1, title2)

    # Feed to plotting function
    #atm.add_atm_plot(P, T, gas_profiles, molec_names, legend=True, title=plot_title)
    atm.plot_double_atm(P, T, gas_profiles, molec_names, title=title,
                        savetag=savetag, xlim=xlim,
                        ylim=ylim, tlim=tlim, seed=seed, legloc=legloc)

    return
#########################################

if __name__ == "__main__":

    make_fig()
