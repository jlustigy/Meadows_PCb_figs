import numpy as np
import coronagraph as cg
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

__all__ = ["add_coronagraph_axes", "plot_coronagraph"]

def add_coronagraph_axes(lam, sol, rad, alpha, phase=90., diam=30., lammin=0.35, lammax=3.0,
                    ref_lam=0.76, wantsnr=10.0, itime=None, Tput=0.2, saveplot=False, title="30m",
                    iwa_lines=True, ax1=None, ax2=None, Tsys=300., Tdet=50.,
                    ground=False):

    lw = 2.0
    mpl.rcParams['font.size'] = 20.0

    # index of phase
    idx = np.argmin(np.abs(alpha - phase))

    # Select phase
    lam = lam[idx,::-1]
    sol = sol[idx,::-1]
    rad = rad[idx,::-1]

    # Filter NaNs
    mask = np.isfinite(rad) & np.isfinite(sol)
    lam = lam[mask]
    sol = sol[mask]
    rad = rad[mask]

    ################################
    # CORONAGRAPH MODEL
    ################################

    # Planet params
    Phi   = 1.      # phase function at quadrature (already included in SMART run)
    Rp    = 1.074     # Earth radii
    r     = 0.0485     # semi-major axis (AU)
    # Stellar params
    Teff  = 3042.   # Sun-like Teff (K)
    Rs = 0.141      # star radius in solar radii
    # Planetary system params
    d = 1.302      # distance to system (pc)
    Nez  = 1.      # number of exo-zodis

    # Telescope parameters
    #lammin = 0.3
    #lammax = 3.0
    Res    = 70.0
    #diam   = 16.0
    #Tput   = 0.05
    C      = 1e-10
    IWA    = 1.0
    OWA    = 40.0
    #Tsys   = 280.0
    #Tdet   = 50.0
    emis   = 0.9
    De     = 1e-4
    DNHpix = 3.0
    Re     = 0.1
    Dtmax  = 1.0
    X      = 1.5
    qe     = 0.9
    MzV    = 23.0
    MezV   = 22.0

    # set lammax to IWA = 1
    if iwa_lines:
        sep =  r / d # arcsec
        iwa1 = 1.22 * (lam * 1e-6) / diam * 206265 # arcsec
        iwa2 = 2. * (lam * 1e-6) / diam * 206265 # arcsec
        iwa3 = 3. * (lam * 1e-6) / diam * 206265 # arcsec
        idy = np.argmin(np.abs(iwa1 - sep))
        lammax = lam[idy]

    # Calculate hi-resolution reflectivity
    Ahr   = np.pi * rad / sol

    ################################
    # RUN CORONAGRAPH MODEL
    ################################

    # Run coronagraph with default LUVOIR telescope (aka no keyword arguments)
    lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
        cg.count_rates(Ahr, lam, sol, phase, Phi, Rp, Teff, Rs, r, d, Nez,
                       THERMAL = True, wantsnr=wantsnr, GROUND=ground,
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

    # Get reference wavelength
    idz = np.argmin(np.abs(lam - ref_lam))
    Dt = DtSNR[idz]
    if itime is not None:
        Dt = itime

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

    # Create figure
    if (ax1 is None) and (ax2 is None):
        fig = plt.figure(figsize=(10,14))
        gs = gridspec.GridSpec(2,1)
        ax1 = plt.subplot(gs[1])
        ax2 = plt.subplot(gs[0])
        plt.subplots_adjust(wspace=0, hspace=0.07)
    else:
        pass

    xlim =  [lammin-2*dlam[0], lammax+dlam[-1]]

    # Set string for plot text
    if Dt > 2.0:
        timestr = "{:.0f}".format(Dt)+' hours'
    else:
        timestr = "{:.0f}".format(Dt*60)+' mins'
    plot_text = r'Integration time = '+timestr

    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = (np.abs(lam - ref_lam)).argmin()
        ref_SNR = SNR[ireflam]
        plot_text = plot_text + '\n SNR = '+"{:.1f}".format(ref_SNR)+\
            ' at '+"{:.2f}".format(lam[ireflam])+r' $\mu$m'

    # Draw plot
    ax2.plot(lam, Cratio*1e9, lw=2.0, color="purple", alpha=0.7, ls="steps-mid")
    ax2.errorbar(lam, spec*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)

    # Set labels
    ax2.set_ylabel(r"F$_p$/F$_s$ ($\times 10^9$)")
    ax2.set_xlabel("Wavelength [$\mu$m]")
    ax2.text(0.75, 0.97, plot_text,\
         verticalalignment='top', horizontalalignment='center',\
         transform=ax2.transAxes,\
         color='black', fontsize=15, bbox=dict(boxstyle="square", fc="w", ec="k"))

    # Create vertical lines for IWA cutoffs
    if iwa_lines:
        sep =  r / d # arcsec
        iwa1 = 1.22 * (lam * 1e-6) / diam * 206265 # arcsec
        iwa2 = 2. * (lam * 1e-6) / diam * 206265 # arcsec
        iwa3 = 3. * (lam * 1e-6) / diam * 206265 # arcsec
        idy1 = np.argmin(np.abs(iwa1 - sep))
        idy2 = np.argmin(np.abs(iwa2 - sep))
        idy3 = np.argmin(np.abs(iwa3 - sep))
        ax2.axvline(x = lam[idy1], color="red", lw=2.0, ls="--")
        ax2.axvline(x = lam[idy2], color="green", lw=2.0, ls="--")
        ax2.axvline(x = lam[idy3], color="blue", lw=2.0, ls="--")
        ax1.axvline(x = lam[idy1], color="red", lw=2.0, ls="--")
        ax1.axvline(x = lam[idy2], color="green", lw=2.0, ls="--")
        ax1.axvline(x = lam[idy3], color="blue", lw=2.0, ls="--")

    # Adjust x,y limits
    ylim = [-10, np.max(Cratio)*1e9+20]
    if ylim is not None: ax2.set_ylim(ylim)
    if xlim is not None: ax2.set_xlim(xlim)

    #####################################

    # Create figure
    #ax = plt.subplot(gs[0])

    # Draw plot
    if itime is None:
        ax1.plot(lam, DtSNR, lw=2.0, color="purple", alpha=0.7, ls="steps-mid")
        ax1.semilogy()
        ax1.set_ylabel(r"Integration Time [hours]")
        # Adjust x,y limits
        ylim = [np.min(DtSNR[np.nonzero(DtSNR)]), 8760.]
        if ylim is not None: ax1.set_ylim(ylim)
        if xlim is not None: ax1.set_xlim(xlim)
    else:
        ax1.plot(lam, SNR, lw=2.0, color="purple", alpha=0.7, ls="steps-mid")
        ax1.semilogy()
        ax1.set_ylabel(r"SNR")
        # Adjust x,y limits
        ylim = [1.0, np.max(SNR[np.nonzero(SNR)])]
        if ylim is not None: ax1.set_ylim(ylim)
        if xlim is not None: ax1.set_xlim(xlim)

    #ax.set_title(title)

    return ax1, ax2

def plot_coronagraph(alpha, output, wantsnr=20.0, itime=None, savetag="fig", plotdir="../../figures/"):

    lammin = 0.3
    lammax = 1.5
    Tput = 0.05
    amin = 0.0
    amax = 150.0
    lw = 2.0

    i=0

    lam = output[i][0][:,:]
    sol = output[i][1][:,:]
    rad = output[i][2][:,:]
    alb = output[i][3][:,:]

    fig = plt.figure(figsize=(20,10))
    gs = gridspec.GridSpec(2,3)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[3])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[4])
    ax4 = plt.subplot(gs[2])
    ax5 = plt.subplot(gs[5])
    plt.subplots_adjust(wspace=0.07, hspace=0.07)


    # Plot luvoir
    ax2,ax3 = add_coronagraph_axes(lam, sol, rad, alpha,
                               Tput=Tput, diam=16.0,
                               title="LUVOIR 16m",
                               ax1=ax2, ax2=ax3, Tsys=200.0, Tdet=50.0,
                               ground=False, wantsnr=wantsnr, itime=itime)
    ax2.set_title("LUVOIR 16m")

    # Plot HabEx
    ax0, ax1 = add_coronagraph_axes(lam, sol, rad, alpha,
                                Tput=Tput, diam=6.5,
                                title="HabEx 6.5m",
                                ax1=ax0, ax2=ax1, Tsys=200.0, Tdet=50.0,
                                ground=False, wantsnr=wantsnr, itime=itime)
    ax0.set_title("HabEx 6.5m")

    # Plot 30m Ground w/IWA lines
    ax4, ax5 = add_coronagraph_axes(lam, sol, rad, alpha,
                                       Tput=Tput, title="Ground-Based 30m",
                                       lammin = 0.35, lammax = 1.5,
                                       ax1=ax4, ax2=ax5, iwa_lines=True,
                                       Tsys=269.0, Tdet=50.0, ground=True,
                                       wantsnr=wantsnr, itime=itime
                                       )
    ax4.set_title("Ground-Based 30m")

    # Plot 30m Ground
    #fig = plot_coronagraph_ground(lam, sol, rad, Tput=0.05, title="Ground-Based 30m", wantsnr=20, lammin=0.3, lammax=1.5,
    #                              iwa_lines=False)

    # Save plot
    fig.tight_layout()
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".png"), bbox_inches='tight')
    print "Saved:", savetag
