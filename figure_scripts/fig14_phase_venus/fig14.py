"""
Figure 14: Clear-sky Venus-like Phase Curves

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    import numpy as np
    import os
    import sys
    sys.path.append("../../")
    from utils import fig_params
    import utils.phasecurves as pcs

    # Params specific to this plot
    savetag = "fig14"
    lammin = 6.5
    lammax = 18.
    R = 4
    iout = 1

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs/")
    alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])

    # Read-in all disk integrated spectra
    typedir = "phaseout_thermal_noTcontrast/"
    output1, flist = pcs.open_phase_dir(alpha, planetdir, typedir)

    typedir = "phaseout_nothermal/"
    output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)

    typedir = "phaseout_thermal_20/"
    output3, tmp = pcs.open_phase_dir(alpha, planetdir, typedir)

    # Make plot
    pcs.plot_binned_phasecurves(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag)

    return
#########################################

if __name__ == "__main__":

    make_fig()
