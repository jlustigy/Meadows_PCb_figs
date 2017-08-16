"""
Figure 11: All phasecurve plots combined

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
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams['font.size'] = 30.0

    savetag = "fig11_all_phasecurves"
    alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])

    # Set up plot
    N = 4
    hspace = 0.08
    h_percent = 0.55
    w_percent = 0.3
    fig, ax = plt.subplots(N, 1, figsize=(14,(8+hspace)*N))
    fig.subplots_adjust(hspace=hspace)
    axi_list = []

    label_font_size = 35
    tick_font_size = 25
    leg_font_size = 24

    # Fake loop over simulations

    ################################## EARTH ##################################
    i = 0
    bbox = ax[i].get_position()
    left, bottom, width, height = [bbox.x0, bbox.y0+bbox.height*(1.-h_percent), w_percent, bbox.height*h_percent]
    axi = fig.add_axes([left, bottom, width, height])
    axi_list.append(axi)
    ax0 = (ax[i], axi)

    # Params specific to this plot
    title = "Earth-like"
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
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs1_earth/")

    typedir = key+"/Tcontrast_max_idl/"; iout = 0
    output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    typedir = key+"/Tcontrast_20K_idl/"
    output3, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    # Dummy for backwards compatibility
    output1 = output2

    pcs.plot_binned_phasecurves_miri(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag,
                                legloc=legloc, moleloc=None, title=title,
                                saveplot=False, ax0=ax0)

    ax[i].set_xlabel("")

    ################################## ARCHEAN ##################################
    i = 1
    bbox = ax[i].get_position()
    left, bottom, width, height = [bbox.x0, bbox.y0+bbox.height*(1.-h_percent), 0.3, bbox.height*h_percent]
    axi = fig.add_axes([left, bottom, width, height])
    axi_list.append(axi)
    ax0 = (ax[i], axi)

    # Params specific to this plot
    title = "Archean Earth-like"
    lammin = 6.5
    lammax = 18.
    R = 4

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs2_archean/")

    # Read-in all disk integrated spectra
    tpe = "haze"
    key = "combined"
    iout = 0

    typedir = tpe+"_"+key+"/Tcontrast_max_idl/"
    output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    typedir = tpe+"_"+key+"/Tcontrast_20K_idl/"
    output3, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    output1 = output2

    pcs.plot_binned_phasecurves_miri(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag,
                                legloc=legloc, moleloc=None, title=title,
                                saveplot=False, ax0=ax0)

    ax[i].set_xlabel("")

    ################################## HIGH O2 ##################################
    i = 2
    bbox = ax[i].get_position()
    left, bottom, width, height = [bbox.x0, bbox.y0+bbox.height*(1.-h_percent), 0.3, bbox.height*h_percent]
    axi = fig.add_axes([left, bottom, width, height])
    axi_list.append(axi)
    ax0 = (ax[i], axi)

    # Params specific to this plot
    title = r"10 bar O$_2$ (Desiccated)"
    lammin = 6.5
    lammax = 26.3
    R = 3
    iout = 1
    legloc = (0.4, 0.72)

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs3_highO2/")

    # Read-in all disk integrated spectra
    key = "10bar_dry"
    iout = 0

    typedir = key+"/max/"
    output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    typedir = key+"/20K/"
    output3, flist = pcs.open_phase_dir(alpha, planetdir, typedir)
    print flist[iout]

    output1 = output2

    pcs.plot_binned_phasecurves_miri(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag,
                                legloc=legloc, moleloc=None, title=title,
                                saveplot=False, ax0=ax0)

    ax[i].set_xlabel("")

    ################################## VENUS ##################################
    i = 3
    bbox = ax[i].get_position()
    left, bottom, width, height = [bbox.x0, bbox.y0+bbox.height*(1.-h_percent), 0.3, bbox.height*h_percent]
    axi = fig.add_axes([left, bottom, width, height])
    axi_list.append(axi)
    ax0 = (ax[i], axi)

    # Params specific to this plot
    title = "90 bar Venus-like"
    lammin = 6.5
    lammax = 18.
    R = 4
    iout = 0
    legloc = (0.4, 0.72)

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs4_venus/")

    # Read-in all disk integrated spectra
    typedir = "Tcontrast_max_idl/"
    output2, flist = pcs.open_phase_dir(alpha, planetdir, typedir)

    typedir = "Tcontrast_20K_idl/"
    output3, tmp = pcs.open_phase_dir(alpha, planetdir, typedir)

    output1 = output2

    pcs.plot_binned_phasecurves_miri(alpha, output1, output2, output3,
                                R=R, lammin=lammin, lammax=lammax,
                                iout=iout, savetag=savetag,
                                legloc=legloc, moleloc=None, title=title,
                                saveplot=False, ax0=ax0)

    ####################################################################
    # PLOT
    ####################################################################

    for i in range(N):
        title = ax[i].get_title()
        ax[i].set_title("")
        bbox = ax[i].get_position()
        ax[i].text(.40, 0.97, title, ha="left", va="top", transform=ax[i].transAxes)
        plt.setp(ax[i].get_xticklabels(), fontsize=tick_font_size, rotation=0)
        plt.setp(ax[i].get_yticklabels(), fontsize=tick_font_size, rotation=0)
        plt.setp(axi_list[i].get_xticklabels(), fontsize=tick_font_size, rotation=0)
        plt.setp(axi_list[i].get_yticklabels(), fontsize=tick_font_size, rotation=0)
        if i>0:
            ax[i].legend_.remove()
        else:
            ax[i].legend(loc="lower right", fontsize=leg_font_size, ncol=1, columnspacing=0,
                         frameon = False)
            #leg=ax[i].legend(loc=(0.4, 0.72), fontsize=16, ncol=2)
            #leg.get_frame().set_alpha(0.0)

    plt.rcParams.update({'axes.titlesize': label_font_size})
    plt.rcParams.update({'axes.labelsize': label_font_size})


    #import pdb; pdb.set_trace()
    #leg=ax[0].legend(loc=0, fontsize=16, ncol=2)
    #leg.get_frame().set_alpha(0.0)
    #lns = ax[0].get_lines()
    #labs = [l.get_label() for l in lns]
    #ax[0].legend(lns, labs, loc=0)

    fig.savefig("../../figures/"+savetag+".pdf", bbox_inches="tight")
    fig.savefig("../../figures/"+savetag+".eps", bbox_inches="tight")


    return
#########################################

if __name__ == "__main__":

    make_fig()
