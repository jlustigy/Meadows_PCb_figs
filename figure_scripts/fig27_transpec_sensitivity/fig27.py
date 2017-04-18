"""
Figure 27:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig27"
    fname3 = "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat.trnst"
    fname1 = "profile_Earth_proxb_.pt_hitran2012_o4_noh2o_noh2co_187Kstrat.trnst"
    fname2 = "profile_Earth_proxb_.pt_hitran2012_o4_noch4_noh2co_187Kstrat.trnst"
    title3 = r"Earth-like: H$_2$O and CH$_4$ Sensitivity"
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
    plot_kwargs = [
                    {"color" : "C0", "label" : "H$_2$O Contribution Removed", "alpha" : 0.7},
                    {"color" : "orange", "label" : "CH$_4$ Contribution Removed", "alpha" : 0.7},
                    {"color" : "black", "label" : "Earth-Like: Full Spectrum", "alpha" : 0.7}
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
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/", f) for f in [fname1, fname2, fname3]]
    title = [title1, title2, title3]
    # Make reflectance plot
    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=["","",""], labels=None, moleloc=moleloc, forced_single=True,plot_kwargs=plot_kwargs, legend=True)

    return
#########################################

if __name__ == "__main__":

    make_fig()
