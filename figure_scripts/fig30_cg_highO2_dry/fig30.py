"""
Figure 30: Dry high O2 coronagraph plot

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
    import utils.coronagraphy as pcg

    # Params specific to this plot
    savetag = "fig30_FpFs_1.2hr"
    ytype = "FpFs"
    itime = 1.2    # Exposure time
    wantsnr = 20.0  # Desired signal-to-noise ratio

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs/")
    alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])

    # Read-in all disk integrated spectra
    typedir = ""
    output1, flist = pcs.open_phase_dir(alpha, planetdir, typedir)

    # Make plot
    pcg.plot_coronagraph(alpha, output1, savetag=savetag, itime=itime, wantsnr=wantsnr, ytype=ytype)

    return
#########################################

if __name__ == "__main__":

    make_fig()
