import os
import numpy as np
from scipy import integrate
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors

import coronagraph as cg

def colorize(vector,cmap='plasma', vmin=None, vmax=None):
    """Convert a vector to RGBA colors.

    Parameters
    ----------
    vector : array
        Array of values to be represented by relative colors
    cmap : str (optional)
        Matplotlib Colormap name
    vmin : float (optional)
        Minimum value for color normalization. Defaults to np.min(vector)
    vmax : float (optional)
        Maximum value for color normalization. Defaults to np.max(vector)

    Returns
    -------
    vcolors : np.ndarray
        Array of RGBA colors
    scalarmap : matplotlib.cm.ScalarMappable
        ScalerMap to convert values to colors
    cNorm : matplotlib.colors.Normalize
        Color normalization
    """

    if vmin is None: vmin = np.min(vector)
    if vmax is None: vmax = np.max(vector)

    cm = plt.get_cmap(cmap)
    cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
    scalarmap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    vcolors = scalarmap.to_rgba(vector)

    return vcolors,scalarmap,cNorm

def read_one_phasecurve(pfile):
    """
    Reads-in one smart_phasecurve output txt file.

    Parameters
    ----------
    pfile : str
        Filename of smart_phasecurve output file

    Returns
    -------
    lam : 1d array
        Wavelength grid [microns]
    solspec : 1d array
        Stellar flux spectrum [W/m^2/um]
    rad : 1d array
        Disk-integrated planetary radiance [W/m^2/um/sr]
    alb : 1d array
        Apparent albedo
    """
    try:
        data = np.genfromtxt(pfile, skip_header=4)
        lam = data[:,0]
        solspec = data[:,1]
        rad = data[:,2]
        alb = data[:,3]
        return lam, solspec, rad, alb
    except IOError:
        print "Failed to read file. Check file form and directory."
        return None, None, None, None

def read_phasecurves(alpha, form, fdir):
    """
    Reads-in an ensemble of smart_phasecurve output txt files across a grid
    in phase angle and returns 2d arrays of dimension: number of phase angles by
    number of wavelength points.

    Parameters
    ----------
    alpha : numpy array
        Phase angles corresponding to smart_phasecurve output [degrees]
    form : str
        Unique filename tag used for smart run (ex: form="earth_icrccm")
    dir : str
        Directory containing smart_phasecurve files (ex: dir="smart_output/")

    Returns
    -------
    lam : 2d array
        Wavelength grids [microns]
    solspec : 2d array
        Stellar flux spectra [W/m^2/um]
    rad : 2d array
        Disk-integrated planetary radiances [W/m^2/um/sr]
    alb : 2d array
        Apparent albedos
    """
    for i in range(len(alpha)):
        a=str(alpha[i])
        if np.sign(alpha[i]) < 0:
            a = a.ljust(8, '0')
        else:
            a = a.ljust(7, '0')
        pfile = fdir+"phasecurve_"+form+"_"+a+".txt"
        if i==0:
            lam, solspec, rad, alb = read_one_phasecurve(pfile)
        else:
            tlam, tsolspec, trad, talb = read_one_phasecurve(pfile)
            lam = np.vstack([lam, tlam])
            solspec = np.vstack([solspec, tsolspec])
            rad = np.vstack([rad, trad])
            alb = np.vstack([alb, talb])
    return lam, solspec, rad, alb

def plot_disk_integrated_spec(lam, sol, rad, alb, alpha,
                              lammin=None, lammax=None,
                              amin=None, amax=None,
                              lw=2.0, title=""):
    """
    Summary plot for an ensemble of smart_phasecurve runs at different phases.
    Note that 4 of the 5 arguments are exact outputs from read_phasecurves(),
    while 'alpha' is an argument of read_phasecurves().

    Parameters
    ----------
    lam : array
        Wavelength [microns]
    sol : array
        Stellar flux [W/m^2]
    rad : array
        Disk-integrated planet radiance [W/m^2/um/sr]
    alb : array
        Apparent albedo
    alpha : array
        Phase angles

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure object
    """

    # Define min & max wavelengths if not user specified
    if lammin is None:
        lammin = np.min(lam[0])
    if lammax is None:
        lammax = np.max(lam[0])
    # Define min & max phase angles if not user specified
    if amin is None:
        amin = np.min(alpha)
    if amax is None:
        amax = np.max(alpha)

    # Select wavelength and phase range
    mask = (lam[0,:] >= lammin) & (lam[0,:] <= lammax)
    amask = (alpha >= amin) & (alpha <= amax)

    # Create color grid
    color = colorize(alpha[amask])[0]

    # Create figure
    fig = plt.figure(figsize=(12,16))
    gs = gridspec.GridSpec(3,1)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax0.set_xticklabels([])
    ax1.set_xticklabels([])
    plt.subplots_adjust(wspace=0, hspace=0.1)

    # Set axis labels
    ax0.set_ylabel(r"Stellar Spectrum [W m$^{-2}$ $\mu$m$^{-1}$]")
    ax1.set_ylabel(r"Planet Radiance [W m$^{-2}$ $\mu$m$^{-1}$ sr$^{-1}$]")
    ax2.set_ylabel(r"Apparent Albedo")
    ax2.set_xlabel(r"Wavelength [$\mu$m]")
    ax0.set_title(title)

    # Plot solar spectrum
    ax0.plot(lam[0,mask], sol[0,mask], lw=lw, color="k")

    # Loop over phase angles, plotting each disk-integrated radiance and albedo
    for i in range(len(alpha[amask])):
        ax1.plot(lam[amask,:][i,mask], rad[amask,:][i,mask], lw=lw, color=color[i], label=r"$\alpha = %i ^{\circ}$" %alpha[amask][i])
        ax2.plot(lam[amask,:][i,mask], alb[amask,:][i,mask], lw=lw, color=color[i])

    # Set legend
    leg=ax1.legend(loc=0, fontsize=14, ncol=2)
    leg.get_frame().set_alpha(0.0)

    return fig

def plot_alb_phase(lam, alb, alpha, ax=None,
                   lammin=None, lammax=None,
                   amin=None, amax=None,
                   lw=2.0, title=""):
    """
    Plots Planet apparent albedo vs wavelength at many phases.
    """

    # Define min & max wavelengths if not user specified
    if lammin is None:
        lammin = np.min(lam[0])
    if lammax is None:
        lammax = np.max(lam[0])
    # Define min & max phase angles if not user specified
    if amin is None:
        amin = np.min(alpha)
    if amax is None:
        amax = np.max(alpha)

    # Select wavelength and phase range
    mask = (lam[0,:] >= lammin) & (lam[0,:] <= lammax)
    amask = (alpha >= amin) & (alpha <= amax)

    # Create color grid
    color = colorize(alpha[amask])[0]

    # Create figure, set axis
    if ax is None:
        fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(1,1)
        ax1 = plt.subplot(gs[0])
    else:
        ax1 = ax

    # Set axis labels
    ax1.set_ylabel(r"Apparent Albedo")
    ax1.set_xlabel(r"Wavelength [$\mu$m]")
    ax1.set_title(title)

    # Loop over phase angles, plotting each disk-integrated radiance and albedo
    for i in range(len(alpha[amask])):
        ax1.plot(lam[amask,:][i,mask], alb[amask,:][i,mask], lw=lw, color=color[i], label=r"$\alpha = %i ^{\circ}$" %alpha[amask][i])

    # Set legend
    leg=ax1.legend(loc=0, fontsize=14, ncol=2)
    leg.get_frame().set_alpha(0.0)

    if ax is None:
        return fig
    else:
        return

def plot_rad_phase(lam, rad, alpha, ax=None,
                   lammin=None, lammax=None,
                   amin=None, amax=None,
                   lw=2.0, title="", xlog=False, ylog=False):
    """
    Plots Planet radiance vs wavelength at many phases.
    """

    # Define min & max wavelengths if not user specified
    if lammin is None:
        lammin = np.min(lam[0])
    if lammax is None:
        lammax = np.max(lam[0])
    # Define min & max phase angles if not user specified
    if amin is None:
        amin = np.min(alpha)
    if amax is None:
        amax = np.max(alpha)

    # Select wavelength and phase range
    mask = (lam[0,:] >= lammin) & (lam[0,:] <= lammax)
    amask = (alpha >= amin) & (alpha <= amax)

    # Create color grid
    color = colorize(alpha[amask])[0]

    # Create figure, set axis
    if ax is None:
        fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(1,1)
        ax1 = plt.subplot(gs[0])
    else:
        ax1 = ax

    # Set axis labels
    ax1.set_ylabel(r"Planet Radiance [W m$^{-2}$ $\mu$m$^{-1}$ sr$^{-1}$]")
    ax1.set_xlabel(r"Wavelength [$\mu$m]")
    ax1.set_title(title)

    # Loop over phase angles, plotting each disk-integrated radiance and albedo
    for i in range(len(alpha[amask])):
        ax1.plot(lam[amask,:][i,mask], rad[amask,:][i,mask], lw=lw, color=color[i], label=r"$\alpha = %i ^{\circ}$" %alpha[amask][i])

    # Set legend
    leg=ax1.legend(loc=0, fontsize=14, ncol=2)
    leg.get_frame().set_alpha(0.0)

    # Set additional params
    ax1.set_xlim([lammin, lammax])
    if xlog: ax1.semilogx()
    if ylog: ax1.semilogy()

    if ax is None:
        return fig
    else:
        return

def plot_phasecurve(alpha, wl, form="", fdir="", title=""):
    """
    Plots planet radiance at various wavelengths as a function of phase.
    """

    # Read disk integrated phases
    lam, solspec, rad, alb = read_phasecurves(alpha, form, fdir)

    # Find lam indicies nearest wl grid points
    iwl = np.array([find_nearest(lam[0,:], wlx) for wlx in wl])

    # Get colors for each wavelength point
    color = colorize(np.arange(len(wl)), cmap="viridis")[0]

    # Create figure
    fig = plt.figure(figsize=(12,10))
    gs = gridspec.GridSpec(1,1)
    ax0 = plt.subplot(gs[0])

    # Set labels
    ax0.set_xlabel(r"Phase Angle [$^{\circ}$]")
    ax0.set_ylabel(r"Planet Radiance [W m$^{-2}$ $\mu$m$^{-1}$ sr$^{-1}$]")
    ax0.set_title(title)

    # Plot phasecurves for each band
    for i in range(len(wl)):
        ax0.plot(alpha, rad[:,iwl[i]], "o", lw=2.0, ls="-", color=color[i], label=r"$\lambda = %.2f \ \mu$m" %lam[0,iwl[i]])

    # Make legend
    leg=ax0.legend(loc=0, fontsize=14, ncol=2)
    leg.get_frame().set_alpha(0.0)

    # Set additional parameters
    ax0.set_xlim([np.min(alpha)-5, np.max(alpha)+5])

    return fig

def find_nearest(array,value):
    """Finds index of array nearest to the value
    """
    idx = (np.abs(array-value)).argmin()
    return idx

# Define planck function
def planck(temp, wav):
    """
    Planck bb function.
    """
    h = 6.62607e-34       # Planck constant (J * s)
    c = 2.998e8           # Speed of light (m / s)
    k = 1.3807e-23        # Boltzmann constant (erg / K)
    wav = wav * 1e-6
    return (2. * h * c**2) / (wav**5) / (np.exp(h * c / (wav * k * temp)) - 1.0)

def fix_stellar_flux(lam, sol):
    """
    Extrapolates stellar spectrum with bb beyond observed range.
    """
    # Define some astro quantities
    R_sun  = 6.955e5
    au     = 1.496e8

    # Define some sys params
    T      = 3042.0
    R_star = 0.14

    # Calculate BB intensity
    B = planck(T, lam)

    # Convert to flux at 1AU, scaled by "necessary factor"
    flux_star = 0.0051 *np.pi * np.power(R_star * R_sun, 2.0) * B / (4. * np.pi * np.power(1 * au, 2.))

    # Construct new flux, using original flux where valid
    validran = (lam < 5.5)
    sol_fix = np.hstack([flux_star[~validran], sol[validran]])

    return sol_fix

def open_phase_dir(alpha, planetdir, typedir):
    """
    Opens a grid of phasecurve files.

    Parameters
    ----------
    alpha : ndarray
        Planetary phase grid
    planetdir : str
        Directory for a particular planet
    typedir : str
        Directory containing disk-integrated output files

    Returns
    -------
    list : list
        list of read-in phasecurves
    flist : list
        names of unique simulation names
    """

    # Construct path to data files
    fdir = os.path.join(planetdir, typedir)

    # Get all file names in directory
    listdir = os.listdir(fdir)

    # Filter out hidden files: .*
    igood = []
    for i in range(len(listdir)):
        if not listdir[i].startswith("."):
            igood.append(i)
    listdir = np.array(listdir)
    listdir = list(listdir[np.array(igood)])


    # Extract list of forms from directory
    forms = np.array([listdir[i][11:-(listdir[i][::-1].find("_") + 1)] for i in range(len(listdir))])

    # isolate unique forms
    flist = []
    for i in range(len(forms)):
        # if form is new AND it isn't a hidden file
        if (forms[i] not in flist) and (not forms[i].startswith(".")):
            # append new form to list
            flist.append(forms[i])

    return [read_phasecurves(alpha, form, fdir) for form in flist], flist

def compute_lam(lammin, lammax, Res):
    """
    Construct wavelength grid (from coronagraph model).
    """
    # Set wavelength grid
    lam  = lammin #in [um]
    Nlam = 1
    while (lam < lammax):
        lam  = lam + lam/Res
        Nlam = Nlam +1
    lam    = np.zeros(Nlam)
    lam[0] = lammin
    for j in range(1,Nlam):
        lam[j] = lam[j-1] + lam[j-1]/Res
    Nlam = len(lam)
    dlam = np.zeros(Nlam) #grid widths (um)
    # Set wavelength widths
    for j in range(1,Nlam-1):
        dlam[j] = 0.5*(lam[j+1]+lam[j]) - 0.5*(lam[j-1]+lam[j])
    #widths at edges are same as neighbor
    dlam[0] = dlam[1]
    dlam[Nlam-1] = dlam[Nlam-2]
    return lam, dlam

def cphot(F, lam, A, T, R):
    """
    Computes photon count rates at each wavelength

    Parameters
    ----------
    F : ndarray
        Object flux density grid [W/m**2/um]
    lam : ndarray
        Corresponding wavelength grid [um]
    A : float
        Telescope aperture [m**2]
    T : float
        Telescope throughput
    R : float
        Telescope resolving power (lam/dlam)
    """
    hc = 1.986446e-25  # h*c (kg*m**3/s**2)
    return F * A * T * lam**2 / (hc * R) * 1e-6

def F_photons(Flam, lam):
    """
    Photon flux given spectral flux density and wavelength
    """
    hc = 1.986446e-25  # h*c (kg*m**3/s**2)
    return Flam * lam**2 * 1e-6 / hc

def phase_phots(lam_in, sol_in, rad_in, alb_in, alpha,
              lammin = 5.5,
              lammax = 19.,
              amax = 180.,
              amin = 0.,
              R = 20.):
    """
    Performs a handful of potentially useful phasecurve calculations.

    Parameters
    ----------
    lam_in : ndarray
        Wavelength grid [um]
    sol_in : ndarray
        Solar flux grid [W/m**2/um]
    rad_in : ndarray
        Planetary radiance grid [W/m**2/um/sr]
    alb_in : ndarray
        Planetary albedo grid
    alpha : ndarray
        Planetary phase grid [degrees]
    lammin : float (optional)
        Minimum wavelength [um]
    lammax : float (optional)
        Maximum wavelength [um]
    amin : float (optional)
        Minimum phase angle [deg]
    amax : float (optional)
        Maximum phase angle [deg]
    R : float (optional)
        Telescope resolving power (lam/dlam)

    Returns
    -------
    lamlo : ndarray
        Low-res wavelength grid
    dlamlo : ndarray
        Low-res wavelength element widths
    Nsarr : ndarray
        Number of stellar photons (depreciated)
    Nparr : ndarray
        Number of planetary photons (depreciated)
    Fsp : ndarray
        Stellar photon flux
    Fpp : ndarray
        Planetary photon flux
    Fs : ndarray
        Stellar flux
    Fp : ndarray
        Planetary flux
    [lam, Fphi, Fshi] : list of ndarrays
        High-res planet and star flux at earth
    """

    # Create wavelength grid
    lamlo, dlamlo = compute_lam(lammin, lammax, R)

    # Set sys params
    a = 0.0485  # AU
    d = 1.302   # pc
    Rp = 1.074  # 6849 km
    Rearth = 6378. # km/Rearth
    au = 1.496e8   # km/AU
    pc = 3.086e13  # km/pc

    # Create lists to hold data
    Fp = np.zeros([len(alpha), len(lamlo)])
    Fs = np.zeros([len(alpha), len(lamlo)])
    Fpp = np.zeros([len(alpha), len(lamlo)])
    Fsp = np.zeros([len(alpha), len(lamlo)])
    Nsarr = np.zeros([len(alpha), len(lamlo)])
    Nparr = np.zeros([len(alpha), len(lamlo)])
    cplan = np.zeros([len(alpha), len(lamlo)])
    cstar = np.zeros([len(alpha), len(lamlo)])
    Fphi = np.zeros([len(alpha), len(lam_in[0,:])])
    Fshi = np.zeros([len(alpha), len(lam_in[0,:])])

    for j in range(len(alpha)):

        # Select phase
        phase = alpha[j]
        lam = lam_in[j,:]
        sol = sol_in[j,:]
        rad = rad_in[j,:]
        alb = alb_in[j,:]

        # Use stellar BB beyond 5.5 microns
        sol = fix_stellar_flux(lam, sol)

        # Filter NaNs
        #nanmask = np.isfinite(rad) & np.isfinite(sol)
        #lam = lam[nanmask]
        #sol = sol[nanmask]
        #rad = rad[nanmask]
        #alb = alb[nanmask]

        # Energy Flux
        Ftoa = rad * np.pi # W/m^2/um
        Fp_earth = Ftoa * (Rp * Rearth)**2. / ((d * pc)**2.)
        Fs_earth = sol * ((a * au)/(d * pc))**2.
        Fps_earth = Fp_earth + Fs_earth # W/m^2/um

        Fphi[j,:] = Fp_earth
        Fshi[j,:] = Fs_earth

        # photon flux
        Fp_phot = F_photons(Fp_earth, lam)
        Fs_phot = F_photons(Fs_earth, lam)
        Fps_phot = Fp_phot + Fs_phot

        # Exposure time
        period = 11.186 # days
        Nobs = len(alpha)
        texp = period * 24. * 3600. / Nobs # seconds

        # Total photons
        A = 25. # Collecting Area [m**2]
        T = 1.0 # Throughput
        Np = Fp_phot * A * texp / R
        Ns = Fs_phot * A * texp / R

        # Integrate photons at high resolution
        mask = (lam > (lammin - 0.5*lammin/R)) & (lam < (lammax + 0.5*lammax/R))
        iNp_hr = integrate.trapz(Np[mask][::-1], x=lam[mask][::-1])
        iNs_hr = integrate.trapz(Ns[mask][::-1], x=lam[mask][::-1])


        # Degrade spectrum
        Fp[j,:] = cg.degrade_spec(Fp_earth, lam, lamlo, dlam=dlamlo)
        Fs[j,:] = cg.degrade_spec(Fs_earth, lam, lamlo, dlam=dlamlo)
        Fpp[j,:] = cg.degrade_spec(Fp_phot, lam, lamlo, dlam=dlamlo)
        Fsp[j,:] = cg.degrade_spec(Fs_phot, lam, lamlo, dlam=dlamlo)
        Nparr[j,:] = cg.degrade_spec(Np, lam, lamlo, dlam=dlamlo)
        Nsarr[j,:] = cg.degrade_spec(Ns, lam, lamlo, dlam=dlamlo)

        # Integrate photons at low resolution
        iNp_lr = integrate.trapz(Nparr[j,:], x=lamlo)
        iNs_lr = integrate.trapz(Nsarr[j,:], x=lamlo)

        #print "Ratio of integrated planet photons (hr/lr) =", iNp_hr / iNp_lr
        #print "Ratio of integrated stellar photons (hr/lr) =", iNs_hr / iNs_lr
        #print "-----------------"

        # photon count rate
        T = 1.0 # Throughput
        cplan[j,:] = cphot(Fp[j,:], lamlo, A, T, R)
        cstar[j,:] = cphot(Fs[j,:], lamlo, A, T, R)

    return lamlo, dlamlo, Nsarr, Nparr, Fsp, Fpp, Fs, Fp, [lam, Fphi, Fshi]

def plot_binned_phasecurves(alpha, output1, output2, output3, iout=0, savetag="fig",
                            amin=0.0, amax=180.0, iout20=0, plotdir="../../figures/",
                            R=3, lammin=6.5, lammax=26.3):
    """
    Creates binned phasecurve plots (used in Meadows et al. paper)

    Parameters
    ----------
    alpha : ndarray
        Planetary phase grid [degrees]
    output1 : list
        No T contrast phasecurves read-in from open_phase_dir()
    output2 : list
        No nightside flux phasecurves read-in from open_phase_dir()
    output3 : list
        20K cooler nightside phasecurves read-in from open_phase_dir()
    iout : int (optional)
        index corresponding to which simulation if multiple exist in dir
    iout20 : int (optional)
        index corresponding to which simulation if multiple exist in dir for 20K
    savetag : str
        naming convention for saved figures
    lammin : float (optional)
        Minimum wavelength [um]
    lammax : float (optional)
        Maximum wavelength [um]
    amin : float (optional)
        Minimum phase angle [deg]
    amax : float (optional)
        Maximum phase angle [deg]
    R : float (optional)
        Telescope resolving power (lam/dlam)
    plotdir : str (optional)
        Relative location to save plots
    """

    lamlo1, dlamlo1, Nsarr1, Nparr1, Fsp1, Fpp1, Fs1, Fp1, hires1 = \
    phase_phots(output1[iout][0],output1[iout][1],output1[iout][2],output1[iout][3], \
                alpha, R=R, lammin=lammin, lammax=lammax)

    lamlo2, dlamlo2, Nsarr2, Nparr2, Fsp2, Fpp2, Fs2, Fp2, hires2 = \
    phase_phots(output2[iout][0],output2[iout][1],output2[iout][2],output2[iout][3], \
                alpha, R=R, lammin=lammin, lammax=lammax)

    # No "iout here b/c currently only one case per planet type"
    lamlo3, dlamlo3, Nsarr3, Nparr3, Fsp3, Fpp3, Fs3, Fp3, hires3 = \
    phase_phots(output3[iout20][0],output3[iout20][1],output3[iout20][2],output3[iout20][3], \
                alpha, R=R, lammin=lammin, lammax=lammax)

    lamlo = lamlo1
    dlamlo = dlamlo1

    lamhi1 = hires1[0]
    Fpe1 = hires1[1]
    Fse1 = hires1[2]

    lamhi2 = hires2[0]
    Fpe2 = hires2[1]
    Fse2 = hires2[2]

    lamhi3 = hires3[0]
    Fpe3 = hires3[1]
    Fse3 = hires3[2]

    xaxlammin = lamlo[0] - dlamlo[0] - dlamlo[0]/10
    xaxlammax = lamlo[-1] + dlamlo[-1] + dlamlo[-1]/10

    yaxlammin = lamlo[0] - dlamlo[0]/2
    yaxlammax = lamlo[-1] + dlamlo[-1]/2

    lmask = (lamhi2 > xaxlammin) & (lamhi2 < xaxlammax)
    lmask2 = (lamhi2 > yaxlammin) & (lamhi2 < yaxlammax)
    amask = (alpha >= amin) & (alpha <= amax)

    colors1 = colorize(np.arange(len(lamlo)))[0]
    colors2 = colorize(np.arange(len(alpha[amask])), cmap="viridis")[0]

    fig = plt.figure(figsize=(20,8))
    gs = gridspec.GridSpec(1,2)
    plt.subplots_adjust(wspace=0.05, hspace=0.0)

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    # Set axis labels
    ax0.set_ylabel(r"Planet/Star Contrast")
    ax0.set_xlabel(r"Phase [deg]")
    ax0.semilogy()
    for i in range(len(lamlo)):
        ax0.plot(alpha, Fp1[:,i]/Fs1[:,i], lw=2.0, c=colors1[i], ls="-", label=r"$%.1f\mu$m" % lamlo[i])
        ax0.plot(alpha, Fp2[:,i]/Fs2[:,i], lw=2.0, c=colors1[i], ls="--")
        ax0.plot(alpha, Fp3[:,i]/Fs3[:,i], lw=2.0, c=colors1[i], ls="dotted")
        #ax0.fill_between(alpha, Fp1[:,i]/Fs1[:,i], Fp2[:,i]/Fs2[:,i], color=colors1[i], alpha=0.1)
        ax1.axvspan(lamlo[i] - 0.5*dlamlo[i], lamlo[i] + 0.5*dlamlo[i], color=colors1[i], alpha=0.1)
        ax1.axvline(lamlo[i] - 0.5*dlamlo[i], c=colors1[i], ls="-", lw=1.0)
        ax1.axvline(lamlo[i] + 0.5*dlamlo[i], c=colors1[i], ls="-", lw=1.0)
    leg=ax0.legend(loc=0, fontsize=16, ncol=2)
    leg.get_frame().set_alpha(0.0)

    """
    # Plot all phases
    yarr= Fpe3[amask,:][:,lmask]/Fse3[amask,:][:,lmask]
    ylabel = r"Planet/Star Contrast"
    # Set axis labels
    ax1.set_ylabel(ylabel, rotation=270, labelpad=30)
    ax1.set_xlabel(r"Wavelength [$\mu$m]")
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    for i in range(len(alpha[amask])):
        ax1.plot(lamhi3[lmask], yarr[i,:], lw=2.0, c=colors2[i], label=r"$%.0f ^{\circ}$" % alpha[amask][i])
    ax1.semilogy()
    leg=ax1.legend(loc=0, fontsize=16, ncol=2)
    leg.get_frame().set_alpha(0.0)
    """
    # Plot full phase
    yarr= Fpe2[:,lmask]/Fse2[:,lmask]
    ylabel = r"Planet/Star Contrast"
    # Set axis labels
    ax1.set_ylabel(ylabel, rotation=270, labelpad=30)
    ax1.set_xlabel(r"Wavelength [$\mu$m]")
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    i = find_nearest(alpha, 0.0)
    ax1.plot(lamhi2[lmask], Fpe2[i,lmask]/Fse2[i,lmask], lw=2.0, c="k")
    ax1.semilogy()

    ymin = np.nanmin(Fpe2[i,lmask2]/Fse2[i,lmask2])
    ymax = np.nanmax(Fpe2[i,lmask]/Fse2[i,lmask])

    ax1.set_xlim([xaxlammin, xaxlammax])
    #ymin = np.nanmin(Fpe2[amask,:][:,lmask2]/Fse2[amask,:][:,lmask2])
    #ymax = np.nanmax(Fpe2[amask,:][:,lmask2]/Fse2[amask,:][:,lmask2])
    ax0.set_ylim([ymin, ymax])
    ax1.set_ylim([ymin, ymax])


    # Save plot
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')
    print "Saved:", savetag

    return
