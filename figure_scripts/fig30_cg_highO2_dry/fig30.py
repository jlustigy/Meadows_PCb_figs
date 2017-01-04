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
    savetag = "fig30"

    # More general params
    planetdir = os.path.join(os.path.dirname(__file__),"model_outputs/")
    alpha = np.array([-180, -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])

    # Read-in all disk integrated spectra
    typedir = ""
    output1, flist = pcs.open_phase_dir(alpha, planetdir, typedir)

    # Make plot
    pcg.plot_coronagraph(alpha, output1, savetag=savetag, itime=20., wantsnr=20.0)

    """
    ###
    import coronagraph as cg
    import readsmart as rs
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib import rc
    import pdb
    import sys
    radfile = "profile_o2lb_10bar_dry.pt_filtered_hitran2012_50_100000cm_toa.rad"
    itime = 1.0

    # Planet params
    Phi   = 0.318     # phase function at quadrature (already included in SMART run)
    Rp    = 1.074     # Earth radii
    r     = 0.0485     # semi-major axis (AU)

    # Stellar params
    Teff  = 3042.   # Sun-like Teff (K)
    Rs = 0.141      # star radius in solar radii
    # Planetary system params
    d = 1.302      # distance to system (pc)
    Nez  = 1.      # number of exo-zodis

    # Telescope parameters
    Res    = 70.0
    diam   = 16.0
    lammin = 0.35
    lammax = 2.0
    Tput = 0.05
    C      = 1e-10
    IWA    = 1.0
    OWA    = 40.0
    Tsys   = 200.0
    Tdet   = 50.0
    emis   = 0.9
    De     = 1e-4
    DNHpix = 3.0
    Re     = 0.1
    Dtmax  = 1.0
    X      = 1.5
    qe     = 0.9
    MzV    = 23.0
    MezV   = 22.0

    # Integration time (hours)
    Dt = itime

    # Planet params
    alpha = 90.     # phase angle at quadrature

    # Plot params
    plot = True
    ref_lam = 0.76
    title = ""
    ylim =  None#[-0.1, 0.8]
    xlim =  None
    tag = ""

    # Save params
    savefile = False
    saveplot = True


    ################################
    # READ-IN DATA
    ################################

    # Read-in spectrum file
    fn = os.path.join(os.path.dirname(__file__),radfile)

    # Read-in .rad file
    lamhr, wno, solar_spec, TOA_flux, rad_streams = rs.rad(radfile,getdata=True)

    # Calculate Hi-res reflectivity spectrum
    Ahr = (TOA_flux / solar_spec) #* np.pi / planet.Phi

    solhr = solar_spec

    sep =  r / d # arcsec
    iwa1 = 1.22 * (lamhr * 1e-6) / diam * 206265 # arcsec
    iwa2 = 2. * (lamhr * 1e-6) / diam * 206265 # arcsec
    iwa3 = 3. * (lamhr * 1e-6) / diam * 206265 # arcsec
    idy = np.argmin(np.abs(iwa1 - sep))
    lammax = lamhr[idy]

    ################################
    # RUN CORONAGRAPH MODEL
    ################################

    # Run coronagraph with default LUVOIR telescope (aka no keyword arguments)
    lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
        cg.count_rates(Ahr, lamhr, solhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez,
                       THERMAL = True, wantsnr=20.0, GROUND=False,
                        lammin = lammin,
                        lammax = lammax,
                        Res    = Res   ,
                        diam   = diam  ,
                        Tput   = Tput  ,
                        C      = C     ,
                        IWA    = IWA   ,
                        OWA    = OWA   ,
                        Tsys   = Tsys  ,
                        Tdet   = Tdet  ,
                        emis   = emis  ,
                        De     = De    ,
                        DNHpix = DNHpix,
                        Re     = Re    ,
                        Dtmax  = Dtmax ,
                        X      = X     ,
                        qe     = qe    ,
                        MzV    = MzV   ,
                        MezV   = MezV
                      )

    # Calculate background photon count rates
    cb = (cz + cez + csp + cD + cR + cth)

    # Convert hours to seconds
    Dts = Dt * 3600.

    # Calculate signal-to-noise assuming background subtraction (the "2")
    SNR  = cp*Dts/np.sqrt((cp + 2*cb)*Dts)

    # Calculate 1-sigma errors
    sig= Cratio/SNR

    # Add gaussian noise to flux ratio
    spec = Cratio + np.random.randn(len(Cratio))*sig

    ################################
    # PLOTTING
    ################################

    if plot:

        # Create figure
        fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(1,1)
        ax = plt.subplot(gs[0])

        # Set string for plot text
        if Dt > 2.0:
            timestr = "{:.0f}".format(Dt)+' hours'
        else:
            timestr = "{:.0f}".format(Dt*60)+' mins'
        plot_text = r'Distance = '+"{:.1f}".format(d)+' pc'+\
        '\n Integration time = '+timestr

        # If a reference wavelength is specified then return the SNR at that wl
        # corresponding to the integration time given
        if ref_lam:
            ireflam = (np.abs(lam - ref_lam)).argmin()
            ref_SNR = SNR[ireflam]
            plot_text = plot_text + '\n SNR = '+"{:.2f}".format(ref_SNR)+\
                ' at '+"{:.2f}".format(lam[ireflam])+r' $\mu$m'

        # Draw plot
        ax.plot(lam, Cratio*1e9, lw=2.0, color="purple", alpha=0.7, ls="steps-mid")
        ax.errorbar(lam, spec*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)

        # Set labels
        ax.set_ylabel(r"F$_p$/F$_s$ ($\times 10^9$)")
        ax.set_xlabel("Wavelength [$\mu$m]")
        ax.set_title(title)
        ax.text(0.99, 0.99, plot_text,\
             verticalalignment='top', horizontalalignment='right',\
             transform=ax.transAxes,\
             color='black', fontsize=20)

        # Adjust x,y limits
        if ylim is not None: ax.set_ylim(ylim)
        if xlim is not None: ax.set_xlim(xlim)

        # Save plot if requested
        if saveplot:
            plot_tag = "luvoir_demo.pdf"
            fig.savefig(plot_tag)
            print 'Saved: ' + plot_tag
        else:
            plt.show()
    ###
    """

    return
#########################################

if __name__ == "__main__":

    make_fig()
