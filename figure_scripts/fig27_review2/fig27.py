"""
Figure 27:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig27_new"
    #fname3 = "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat.trnst"
    #fname1 = "profile_Earth_proxb_.pt_hitran2012_o4_noh2o_noh2co_187Kstrat.trnst"
    #fname2 = "profile_Earth_proxb_.pt_hitran2012_o4_noch4_noh2co_187Kstrat.trnst"
    fname3 = "profile_Earth_proxb_.pt_hitran2012_187Kstrat.trnst"
    fname1 = "profile_Earth_proxb_.pt_hitran2012_noh2o_187Kstrat.trnst"
    fname2 = "profile_Earth_proxb_.pt_hitran2012_noch4_187Kstrat.trnst"
    title3 = r"Earth-like"
    title1 = "H$_2$O Contribution Removed"
    title2 = "CH$_4$ Contribution Removed"
    lammin = 0.2
    lammax = 20.
    ylim = [0.0, 100.0]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]
    moleloc = {"CH$_4$" : [(0.9, 31.), (1.16,40.),(1.4, 43.), (2.33,49.), (3.4, 60), (6.6,39.), (7.9,53.)],
                "H$_2$O" : [(0.86, 7.), (5.8, 5.)],
                 }

    moleloc = {
                "O$_2$" : [(0.215, 57.), (0.69, 24.), (0.76, 27.)],
                "O$_3$" : [(0.28, 60.),(0.6, 25.),(4.8,43.),(9.7,39.)],
                "CH$_4$" : [(0.9, 33.), (1.0, 26.), (1.16,40.),(1.4, 43.), (2.33,49.), (3.4, 60), (6.6,39.), (7.9,53.)],
                "H$_2$O" : [(5.8, 34.)],
                "CO$_2$" : [(1.92, 40.), (2.8, 51.), (4.3, 61.),(15.5, 65.)]
            }

    plot_kwargs = [
                    {"color" : "black", "label" : "Full Spectrum", "alpha" : 0.7, "zorder" : 100},
                    {"color" : "C0", "label" : "H$_2$O Contribution Removed", "alpha" : 0.7},
                    {"color" : "orange", "label" : "CH$_4$ Contribution Removed", "alpha" : 0.7},
                           ]

    # import basic dependencies
    import numpy as np
    import os
    import sys

    # import local package
    sys.path.append("../../")
    from utils import fig_params
    import utils.spectra as spc

    # More general params
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/", f) for f in [fname3, fname1, fname2]]
    title = [title1, title2, title3]
    # Make reflectance plot
    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=[title3,"",""], labels=None, moleloc=moleloc, forced_single=True,plot_kwargs=plot_kwargs, legend=True)

    return
#########################################

if __name__ == "__main__":

    make_fig()
