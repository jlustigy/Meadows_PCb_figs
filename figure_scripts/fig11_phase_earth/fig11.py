"""
Figure 11: Earth-like Phase Curves

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
    savetag = "fig11"
    key = "combined"
    lammin = 6.5
    lammax = 26.3
    R = 3
    iout = 1
    legloc = (0.4, 0.72)
    moleloc = {
        "CO$_2$" : [(10.4, 30.0), (15.0, 40.0)],
        "H$_2$O" : [(23.75, 90.0)],
        "O$_3$" : [(9.6, 20.0)]
    }

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs/")
    #alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])
    alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])

    # Read-in all disk integrated spectra
    typedir = "pre_Tcontrast_min/"; iout = 0
    output1, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    #typedir = "phaseout_nothermal/"
    #typedir = "Tcontrast_max_haze/"
    #typedir = "pre_Tcontrast_max/"
    #typedir = "Tcontrast_max_haze_test/"
    typedir = "new/"+key+"/Tcontrast_max_idl/"; iout = 0
    output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    #typedir = "phaseout_thermal_20/"
    #typedir = "Tcontrast_20K_haze/"
    #typedir = "pre_Tcontrast_20K_fix5/"
    #typedir = "Tcontrast_20K_haze_test/"
    typedir = "new/"+key+"/Tcontrast_20K_idl/"
    output3, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    # Read-in all disk integrated spectra
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
                                iout=iout, savetag=savetag+"_orig")
    """

    """
    pcs.plot_binned_phasecurves_new(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag+"_new")
    """

    pcs.plot_binned_phasecurves_miri(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag,
                                legloc=legloc, moleloc=None)

    return
#########################################

if __name__ == "__main__":

    make_fig()
