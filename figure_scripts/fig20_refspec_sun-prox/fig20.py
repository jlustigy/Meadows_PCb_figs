"""
Figure 20:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # import basic dependencies
    import numpy as np
    import os
    import sys

    # import local package
    sys.path.append("../../")
    from utils import fig_params
    import utils.spectra as spc

    # Params specific to this plot
    savetag = "fig20"
    title1 = "Sun/Proxima Comparison"
    plot_kwargs = [
        {"color" : "orange", "label" : "Earth-Sun", "alpha" : 0.7},
        {"color" : "purple", "label" : "Earth-Proxima", "alpha" : 0.7}
    ]
    #title2 = "10 bar O$_2$ (Dessicated)"
    lammin = 0.1
    lammax = 2.5
    ylim = [-0.02, 0.3]
    fname_new1 = "model_outputs/modern_earth_preindustrial_hitran2012_300_100000cm_combined_toa.rad"
    fname_new2 = "model_outputs/profile_Earth_proxb_combined_toa.rad"
    fname_new = [fname_new1, fname_new2]

    fname1 = "modern_earth_preindustrial_hitran2012_300_100000cm_toa.rad"
    fname2 = "modern_earth_preindustrial_hitran2012_300_100000cm_stratocum_toa.rad"
    fname3 = "modern_earth_preindustrial_hitran2012_300_100000cm_cirrus_toa.rad"
    fname4 = "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat_toa.rad"
    fname5 = "profile_Earth_proxb_.pt_stratocum_hitran2012_o4_noh2co_187Kstrat_toa.rad"
    fname6 = "profile_Earth_proxb_.pt_cirrus_hitran2012_o4_noh2co_187Kstrat_toa.rad"
    fname = [
        os.path.join(os.path.dirname(__file__),"model_outputs/", f)
        for f in [fname1, fname2, fname3, fname4, fname5, fname6]
    ]

    moleloc = {
               "O$_3$" : [(0.25, 0.03), (0.5, 0.23)],
               "CH$_4$" : [(0.88, 0.215), (1.0, 0.215), (1.22, 0.20), (1.4, 0.07), (1.65, 0.15), (1.87, 0.07), (2.05,0.1), (2.27,0.12)]
          }

    # Weight each flux
    dic1 = {
        fname[0] : 0.5,
        fname[1] : 0.25,
        fname[2] : 0.25
    }
    dic2 = {
        fname[3] : 0.5,
        fname[4] : 0.25,
        fname[5] : 0.25
    }
    # Weight spectra and create output file
    spc.weight_spectra(dic1, output=fname_new1, ret=False)
    spc.weight_spectra(dic2, output=fname_new2, ret=False)

    # More general params
    title = [title1, ""]

    # Make reflectance plot
    spc.plot_rad(fname_new, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True, moleloc=[moleloc, None])

    return
#########################################

if __name__ == "__main__":

    make_fig()
