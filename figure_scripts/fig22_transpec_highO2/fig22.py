"""
Figure 22:

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    # Params specific to this plot
    savetag = "fig22"
    #fname1 = "profile_o2lb_10bar_h2o.pt_filtered_transit.tran"
    #fname2 = "profile_o2lb_10bar_dry.pt_filtered_transit.tran"
    fname1 = "10bar_O2_wet.pt_filtered_hitran2012_50_100000cm.trnst"
    fname2 = "10bar_O2_dry.pt_filtered_hitran2012_50_100000cm.trnst"
    title1 = "10 bar O$_2$ (Ocean)"
    title2 = "10 bar O$_2$ (Desiccated)"
    lammin = 0.2
    lammax = 20.
    ylim = [0.0, 100.0]
    labels1 = ["O2", "O3", "O4", "CO", "CO2", "H2O"]
    labels2 = ["O2", "O3", "O4", "CO", "CO2"]
    plot_kwargs = [
        {"color" : "purple", "label" : "Ocean", "alpha" : 0.7},
        {"color" : "orange", "label" : "Desiccated", "alpha" : 0.7}
    ]
    m1 = {"O$_2$" : [(0.215, 61.),(0.69, 36.), (0.76, 41.)],
                  "O$_3$" : [(0.28, 87.),(4.8,42.),(9.7,35.)],
                  "O$_4$" : [(0.58,30.5),(0.64,32.),(1.06, 28.), (1.27,38.)],
                  "CO$_2$" : [(1.55, 37.), (2.0, 42.), (2.8, 49.), (4.3, 61.),(15.5, 63.)],
                  "H$_2$O" : [(6.3, 31.)]
                    }
    m2 = {"O$_2$" : [(0.215, 61.), (0.69, 37.5), (0.76, 41.)],
             "O$_3$" : [(0.28, 87.), (0.6, 44.),(4.8,45.),(9.7,50.)],
             "O$_4$" : [(1.06, 28.), (1.27,38.)],
            "CO$_2$" : [(1.55, 37.), (2.0, 42.), (2.8, 49.), (4.3, 62.),(15.5, 63.)]
                    }
    m3 = {"O$_2$" : [(0.215, 61.), (0.69, 37.5), (0.76, 41.)],
             "O$_3$" : [(0.28, 87.), (0.6, 44.),(4.8,45.),(9.7,50.)],
             "O$_4$" : [(0.58,29.5), (0.64,31.), (1.06, 28.), (1.27,38.)],
            "CO$_2$" : [(1.55, 37.), (2.0, 42.), (2.8, 49.), (4.3, 62.),(15.5, 63.)],
            "H$_2$O" : [(6.3, 31.)]
                    }
    moleloc = [m1, m2]

    # import basic dependencies
    import numpy as np
    import os
    import sys

    # import local package
    sys.path.append("../../")
    from utils import fig_params
    import utils.spectra as spc

    # More general params
    fname = [os.path.join(os.path.dirname(__file__),"model_outputs/", f) for f in [fname1, fname2]]
    title = [title1, title2]
    labels = [labels1, labels2]

    # Make reflectance plot
    """
    spc.plot_trans(fname, savetag=savetag+"2", lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, moleloc=moleloc)
    """

    title = [r"10 bar High O$_2$", ""]

    spc.plot_trans(fname, savetag=savetag, lammin=lammin, lammax=lammax, ylim=ylim,
                 title=title, labels=None, forced_single=True,
                 plot_kwargs=plot_kwargs, legend=True, moleloc=m3)

    return
#########################################

if __name__ == "__main__":

    make_fig()
