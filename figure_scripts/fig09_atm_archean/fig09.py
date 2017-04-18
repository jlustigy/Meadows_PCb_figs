"""
Figure 9

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
    savetag = "fig09"
    fname1 = "fig9_left_panel_archean.atm"
    fname2 = "fig9_right_panel_archean.atm"
    title1 = "Archean $1\%$ CH$_4$, $5\%$ CO$_2$"
    title2 = "Archean $1.5\%$ CH$_4$, $5\%$ CO$_2$"
    seed = 7
    xlim = [1e-8, 2]
    ylim = [1e5, 6e-3 ]
    tlim = [145, 305]
    legloc1 = [3e-1, 2e-2, None, None, None, 3e-2, None, None]
    legloc = [legloc1, legloc1]

    # Read in atm files
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname1)
    data1 = np.genfromtxt(atmpath, skip_header=1)
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname2)
    data2 = np.genfromtxt(atmpath, skip_header=1)

    # Parse atm data
    alt1, P1, T1 = data1[:,0], data1[:,1]*1e5, data1[:,2]
    gas_profiles1 = data1[:,3:]
    molec_names1 = ["H$_2$O", "CH$_4$", "O$_2$", "CO$_2$", "C$_2$H$_6$", "CO", "N$_2$"]

    alt2, P2, T2 = data2[:,0], data2[:,1]*1e5, data2[:,2]
    gas_profiles2 = data2[:,3:]
    molec_names2 = ["H$_2$O", "CH$_4$", "O$_2$", "CO$_2$", "C$_2$H$_6$", "CO", "N$_2$"]

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
