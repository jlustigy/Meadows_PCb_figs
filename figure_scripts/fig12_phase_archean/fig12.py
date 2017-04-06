"""
Figure 12: Hazy Archean Phase Curves

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
    savetag = "fig12_idl"
    lammin = 6.5
    lammax = 18.
    R = 4

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs/")
    alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])

    # Read-in all disk integrated spectra
    #typedir = "phaseout_thermal_noTcontrast/"; iout = 2
    #typedir = "Tcontrast_min_haze/"; iout = 3#2
    #typedir = "Tcontrast_min_haze_test/"; iout = 0#2
    tpe = "nohaze"
    key = "combined"
    typedir = tpe+"_"+key+"/Tcontrast_min/"; iout = 0
    output1, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    #typedir = "phaseout_nothermal/"
    #typedir = "Tcontrast_max_haze/"
    typedir = tpe+"_"+key+"/Tcontrast_max_idl/"
    #typedir = "Tcontrast_max_haze_test/"
    output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    #typedir = "phaseout_thermal_20/"
    #typedir = "Tcontrast_20K_haze/"
    typedir = tpe+"_"+key+"/Tcontrast_20K_idl/"
    #typedir = "Tcontrast_20K_haze_test/"
    output3, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    # Make plot
    """
    pcs.plot_binned_phasecurves(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag)
    """

    pcs.plot_binned_phasecurves_miri(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag+"_"+tpe+"_miri_"+key)

    return
#########################################

if __name__ == "__main__":

    make_fig()
