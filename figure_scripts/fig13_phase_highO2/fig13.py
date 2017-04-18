"""
Figure 13: High O2 Phase Curves

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
    savetag = "fig13"
    lammin = 6.5
    lammax = 26.3
    R = 3
    iout = 1
    legloc = (0.4, 0.72)

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs/")
    alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])

    # Read-in all disk integrated spectra
    key = "10bar_dry"
    typedir = key+"/Tcontrast_min/"; iout = 0
    output1, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    #typedir = "phaseout_nothermal/"
    #typedir = "Tcontrast_max_haze/"
    typedir = key+"/max/"
    #typedir = "Tcontrast_max_haze_test/"
    output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    #typedir = "phaseout_thermal_20/"
    #typedir = "Tcontrast_20K_haze/"
    typedir = key+"/20K/"
    #typedir = "Tcontrast_20K_haze_test/"
    output3, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    #typedir = "phaseout_thermal_noTcontrast/"
    #output1, flist = pcs.open_phase_dir(alpha, planetdir, typedir)

    #typedir = "phaseout_nothermal/"
    #output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)

    #typedir = "phaseout_thermal_20/"
    #output3, tmp = pcs.open_phase_dir(alpha, planetdir, typedir)

    # Make plot
    """
    pcs.plot_binned_phasecurves(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag)
    """

    pcs.plot_binned_phasecurves_miri(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag,
                                legloc=legloc)

    return
#########################################

if __name__ == "__main__":

    make_fig()
