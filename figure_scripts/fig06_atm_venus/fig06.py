"""
Figure 6

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
    savetag = "fig06"
    fname1 = "10barCO2_final.atm"
    fname2 = "90barCO2_final.atm"
    title1 = "10 bar"
    title2 = "90 bar"
    seed = 0
    xlim = [1e-8, 2.]
    ylim = [1e7, 1e-2]
    tlim = [155, 700]

    # Read in atm files
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname1)
    data1 = np.genfromtxt(atmpath, skip_header=1)
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname2)
    data2 = np.genfromtxt(atmpath, skip_header=1)

    # Parse atm data
    P1, T1, alt1 = data1[:,0], data1[:,1], data1[:,2]
    gas_profiles1 = data1[:,3:]
    gas_profiles1 = np.hstack([gas_profiles1[:,:5], gas_profiles1[:,7].reshape((len(P1), 1))])
    molec_names1 = ["H$_2$O", "CO$_2$", "SO$_2$", "OCS", "H$_2$SO$_4$", "HCl", "CO", "N$_2$", "H$_2$SO$_4$(orig)"]
    molec_names1 = ["H$_2$O", "CO$_2$", "SO$_2$", "OCS", "H$_2$SO$_4$", "N$_2$"]

    P2, T2, alt2 = data2[:,0], data2[:,1], data2[:,2]
    gas_profiles2 = data2[:,3:]
    gas_profiles2 = np.hstack([gas_profiles2[:,:5], gas_profiles2[:,7].reshape((len(P2), 1))])
    molec_names2 = ["H$_2$O", "CO$_2$", "SO$_2$", "OCS", "H$_2$SO$_4$", "HCl", "OCS", "N$_2$", "H$_2$S", "H$_2$SO$_4$[initial]"]
    molec_names2 = ["H$_2$O", "CO$_2$", "SO$_2$", "OCS", "H$_2$SO$_4$", "N$_2$"]

    # Calculate condensation vapor pressure curve/s
    Pv1, Tv1 = atm.H2SO4_vapor_pandora(P1, T1)
    Pv2, Tv2 = atm.H2SO4_vapor_pandora(P2, T2)

    P = (P1, P2)
    T = (T1, T2)
    gas_profiles = (gas_profiles1, gas_profiles2)
    molec_names = (molec_names1, molec_names2)
    title = (title1, title2)
    cond_curve = ((Pv1, Tv1, "H$_2$SO$_4$"), (Pv2, Tv2, "H$_2$SO$_4$"))

    # Feed to plotting function
    #atm.add_atm_plot(P, T, gas_profiles, molec_names, legend=True, title=plot_title)
    atm.plot_double_atm(P, T, gas_profiles, molec_names, title=title,
                        savetag=savetag, xlim=xlim,
                        ylim=ylim, tlim=tlim, seed=seed, cond_curve=cond_curve)

    return
#########################################

if __name__ == "__main__":

    make_fig()
