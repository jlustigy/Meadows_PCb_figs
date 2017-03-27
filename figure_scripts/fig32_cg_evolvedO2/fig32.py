"""
Figure 32: 10 bar O2/CO2/CO

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
    import utils.spectra as spc

    # Params specific to this plot
    savetag = "fig32"
    frad = "profile_O2_CO2_10bar_prox_hitran2012_50_100000cm_toa.rad"
    ytype = "FpFs"
    itime = 10.    # Exposure time
    wantsnr = 20.0  # Desired signal-to-noise ratio

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs/")
    alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])

    # Read-in all disk integrated spectra
    typedir = ""
    output1, flist = pcs.open_phase_dir(alpha, planetdir, typedir)

    # Make plot
    """
    pcg.plot_coronagraph(alpha, output1, savetag=savetag, itime=itime, wantsnr=wantsnr, ytype=ytype)
    """

    # Plot disk-integrated spectral at 90 deg
    lamhr, sol, rad = pcg.parse_from_phase(output1, alpha, phase=90)
    pcg.plot_coronagraph(lamhr, sol, rad, savetag=savetag+"_v1_disk", itime=itime, wantsnr=wantsnr, ytype=ytype)

    # Using rad file...
    path = os.path.join(os.path.dirname(__file__), frad)
    output = spc.read_rad(path, Numu=4, Nazi=1)
    lamhr, sol, rad = pcg.parse_from_rad(output, phase=90)

    # Plot v1
    pcg.plot_coronagraph(lamhr, sol, rad, savetag=savetag+"_v1_rad", itime=itime, wantsnr=wantsnr, ytype=ytype)
    pcg.plot_coronagraph_mod(lamhr, sol, rad, savetag=savetag+"_v2_rad", itime=itime, wantsnr=wantsnr, ytype=ytype)
    pcg.plot_coronagraph_mod2(lamhr, sol, rad, savetag=savetag+"_v3_rad", itime=itime, wantsnr=wantsnr, ytype=ytype)

    return
#########################################

if __name__ == "__main__":

    make_fig()
