import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors

import pandas as pd
import coronagraph as cg
from scipy.interpolate import interp1d

import sys
sys.path.append("../")
from utils import molecules, fig_params

__all__ = ["read_rad", "plot_rad", "add_trans_plot", "weight_spectra", "plot_trans",
           "read_transit", "plot_transmittance"]

def read_rad(path, Numu=4, Nazi=1):
    """
    General function to open, read, and parse *.rad files from SMART output

    Parameters
    ----------
    path : str
        Path to file with file name
    Numu : int, optional
        Number of upward streams
    Nazi : int, optional
        Number of observer azimuths

    Returns
    -------
    lam : numpy.ndarray
        Wavelength grid [um]
    wno : numpy.ndarray
        Wavenumber grid [1/cm]
    solar : numpy.ndarray
        Stellar flux at planet toa [W/m**2/um]
    toaf : numpy.ndarray
        Top of atmosphere planetary flux [W/m**2/um]
    rads : numpy.ndarray
        Top of atmosphere planetary radiance streams with dimension
        (``Numu*Nazi`` x ``len(lam)``)
    """

    # Number of elements in a row
    Nrow = 4 + Numu*Nazi

    # Convert each line to vector, compose array of vectors
    arrays = np.array([np.array(map(float, line.split())) for line in open(path)])

    # Flatten and reshape into rectangle grid
    arr = np.hstack(arrays).reshape((Nrow, -1), order='F')

    # Parse columns
    lam   = arr[0,:]
    wno   = arr[1,:]
    solar = arr[2,:]
    toaf  = arr[3,:]
    rads  = arr[4:,:]

    return lam, wno, solar, toaf, rads

def add_rad_plot(lam, refl, radius=6850., ax0=None, legend=False,
                 title=None, xlim=None, ylim=None,
                 legloc=None, plot_kwargs={}):
    """
    """

    # Create new figure and axis or use the one provided
    if ax0 is None:
        fig, ax = plt.subplots(figsize=(16,8))
    else:
        ax = ax0

    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel(r"Reflectivity ($\pi I/F$)")

    ax.plot(lam, refl, **plot_kwargs)

    return

def weight_spectra(dic, ret=True, output=None):
    """
    Read multiple rad files and create a weighted mean of the flux and radiances

    Parameters
    ----------
    dic : dict
        Dictionary with rad filenames as the keys and their
        weights as the values

    Returns
    -------
    lam, wno, solar, toaf, rads
    """
    # Create empty lists for quantities to weight
    weighted_fluxes = []
    weighted_rads = []
    # Loop over dictionary of filenames and weights
    for fname, weight in dic.iteritems():
        # Read-in rad file
        lam, wno, solar, toaf, rads = read_rad(fname)
        # weight quantities
        tmp1 = weight * toaf
        tmp2 = weight * rads
        # Append to lists
        weighted_fluxes.append(tmp1)
        weighted_rads.append(tmp2)
    # Convert list to array
    weighted_fluxes = np.array(weighted_fluxes)
    weighted_rads = np.array(weighted_rads)
    # Sum weighted values
    toaf = np.sum(weighted_fluxes, axis=0)
    rads = np.sum(weighted_rads, axis=0)

    # Save output as new rad file
    if output is not None:
        np.savetxt(output, np.vstack([lam, wno, solar, toaf, rads]).T, fmt="%.6e", delimiter="  ")

    # Return values if requested
    if ret:
        return lam, wno, solar, toaf, rads
    else:
        return

def plot_rad(fname, savetag="refl_test", lammin=0.2, lammax=20.0,
                    plot_kwargs={"color" : "black"}, title=None, plotdir="../../figures/",
                    legloc=None, ylim=None, labels=None, forced_single=False,
                    legend=False, Nstr=8, moleloc=None, eps=True):
    """

    Parameters
    ----------
    fname : str or list
        Path + filename of SMART *.rad file
    """

    # Check fname
    if (type(fname) is str) or (type(fname) is unicode):
        # Single rad file provided
        multi = False
        # Make things lists
        fname = [fname]
        title = [title]
        labels = [labels]
        plot_kwargs = [plot_kwargs]
        # Construct figure and axes
        fig, ax = plt.subplots(figsize=(14,6))
    elif type(fname) is list:
        # Multiple rad files provided
        if labels is None: labels = [None for f in fname]
        if type(plot_kwargs) is not list: plot_kwargs = [plot_kwargs for f in fname]
        # Construct figure and axes
        if forced_single:
            multi = False
            fig, ax = plt.subplots(figsize=(14,6))
        else:
            multi = True
            fig, axs = plt.subplots(nrows=len(fname), figsize=(14,(6+1)*len(fname)))
    else:
        print "Invalid fname..."
        return

    # Loop over outputs
    for i in range(len(fname)):
        if multi:
            ax = axs[i]
        f = fname[i]
        label = labels[i]
        # Read in rad file(s)
        lam, wno, solar, toaf, rads = read_rad(f, Numu=int(Nstr/2.))
        mask = (lam > lammin) & (lam < lammax)
        # Calculate reflectivity
        refl = toaf / solar
        # Set labels
        ax.set_xlabel(r"Wavelength [$\mu$m]")
        ax.set_ylabel(r"Reflectivity ($\pi I/F$)")
        ax.set_xlim([lammin, lammax])
        ax.text(0.5, 0.98, title[i],\
             verticalalignment='top', horizontalalignment='center',\
             transform=ax.transAxes,\
             color='black')
        if ylim is not None:
            ax.set_ylim([ylim[0], ylim[1]])
        # Actually plot it!
        ax.plot(lam[mask], refl[mask], **plot_kwargs[i])

        if moleloc is not None:
            if type(moleloc) is list:
                mloc = moleloc[i]
            elif type(moleloc) is dict:
                mloc = moleloc
            else:
                print "Incompatible type for moleloc kwarg"
            if mloc is not None:
                for key, value in mloc.iteritems():
                    #print key, value
                    # get molecule color
                    mcol = molecules.color_from_molecule(key)
                    for im in range(len(value)):
                        # place label
                        ax.text(value[im][0], value[im][1], key, va='center', ha='center',
                             color=mcol, fontsize=16, zorder=10,
                             bbox=dict(boxstyle="square", fc="none", ec="none"))


        # Attempt to label molecules automatically
        if label is not None:
            # Get x and y data
            xdat = ax.lines[0].get_xdata()
            ydat = ax.lines[0].get_ydata()
            # Get x and y axis limits
            xmin, xmax = ax.axes.get_xlim()
            ymin, ymax = ax.axes.get_ylim()

            # degrade spec and put on even grid
            #dlam = (lammax - lammin) / 3000.
            #R = 10
            #xlow, dxlow = cg.noise_routines.construct_lam(lammin, lammax, dlam=dlam)
            #ylow = cg.downbin_spec(ydat, xdat, xlow, dlam=dxlow)
            #ax.plot(xlow, ylow, color="red")

            # Locate rolling max
            #df = pd.DataFrame(refl)
            df = pd.DataFrame(ydat)
            ymaxs = pd.rolling_max(df, int(len(ydat)/15.), min_periods=100)
            #ymaxs = ymaxs[mask].values.reshape(len(xdat))
            ymaxs = ymaxs.values.reshape(len(xdat))
            #ax.plot(xdat, ymaxs, color="blue")

            # Interpolate maxes to lower res
            dlam = (lammax - lammin) / 20.
            xlow, dxlow = cg.noise_routines.construct_lam(lammin, lammax, dlam=dlam)
            f2 = interp1d(xdat, ymaxs, kind='nearest', bounds_error=False, fill_value="extrapolate")
            ynew = f2(xlow)
            #ax.plot(xlow, ynew, color="green")

            # Interpolate back to high res grid
            f2 = interp1d(xlow, ynew, kind='nearest', bounds_error=False, fill_value="extrapolate")
            ycont = f2(xdat)

            # Extract relevant molecules
            molmask = np.array([True if molecules.molecules[j].keys()[0] in label else False
                                for j in range(len(molecules.molecules))])
            mols = molecules.molecules[molmask]
            # Set colors for each molecule
            colors = []
            for i in range(len(mols)):
                if i < 10:
                    ii = i
                else:
                    ii = i - 10
                colors.append("C%i" %ii)
            # Loop over relevant molecules
            for j in range(len(mols)):
                # Extract molecule formula (also dict key)
                key = mols[j].keys()[0]
                # Convert to latex text string
                text = molecules.tex_molecule(key)
                # Loop over bands in molecule
                for k in range(len(mols[j][key])):
                    # is bandcenter in plot window
                    if (mols[j][key][k] > lammin) & (mols[j][key][k] < lammax):
                        # Set x to be location of band
                        x = mols[j][key][k]
                        ix = np.argmin(np.fabs(xdat - x))
                        #y = np.random.uniform(yuse[ix], ymax)
                        y = ycont[ix] + np.random.uniform(0.01, 0.05)
                        #y = ydat[ix]
                        ax.text(x, y, text, va='center', ha='center',
                             color=colors[j], fontsize=12, zorder=100,
                             bbox=dict(boxstyle="square", fc="none", ec="none"))

    # Label
    if legend:
        leg=ax.legend(loc=0, fontsize=16)
        leg.get_frame().set_alpha(0.0)

    # Save figure
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

    if eps:
        fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".eps"), bbox_inches='tight')

    return

def add_trans_plot(lam, solar, toaf, lammin=0.2, lammax=20.0, radius=6850., ax0=None, legend=False,
                 title=None, xlim=None, ylim=None, tlim=None,
                 legloc=None, plot_kwargs={}):
    """
    """

    # Create new figure and axis or use the one provided
    if ax0 is None:
        fig, ax = plt.subplots(figsize=(16,8))
    else:
        ax = ax0

    # Create right y-axis for absorbing radius
    ax1 = ax.twinx()

    #
    mask = (lam > lammin) & (lam < lammax)

    #
    refl = toaf[mask] / solar[mask]

####################################
# TRANSIT
####################################

def read_transit(path, skip_header=0, radius=None):
    """
    Read *.tran or *.trnst output files from SMART.

    Parameters
    ----------
    path : str
        Path/Filename of *.tran file
    radius : tuple, optional
        (planet_radius, stellar_radius) for transit depth calculation [km]

    Returns
    -------
    wl : ndarray
        Wavelength [microns]
    absorbing_radius : ndarray
        Effective absorbing radius of the atmosphere (in km)
        (effective radius - solid body radius)
    tdepth : ndarray
        Transit depth

    Note
    ----
    dF = (Rp/Rs)**2 = -(Rp_solid + column3)**2/Rs**2
    """
    if path.endswith(".tran"):
        data = np.genfromtxt(path, skip_header=0)
        wl = data[:,0]
        flux_ratio = data[:,1]
        absorbing_radius = data[:,2]
        if radius is not None:
            tdepth = (radius[0] + absorbing_radius)**2/radius[1]**2
        else:
            tdepth = None
    elif path.endswith(".trnst"):
        # Read in .trnst file
        data = np.genfromtxt(path, skip_header=2)
        # Split into arrays
        wl = data[:,0]
        wn = data[:,1]
        absorbing_radius = data[:,2]
        tdepth = data[:,3]
    else:
        print "Invalid file..."
        return

    return wl, absorbing_radius, tdepth

def tdepth_from_zeff(zeff, Rp, Rs):
    return (Rp + zeff)**2/Rs**2


def plot_trans(fname, savetag="refl_test", lammin=0.2, lammax=20.0, radius=(6850., 98093.7),
                    plot_kwargs={"color" : "black"}, title=None, plotdir="../../figures/",
                    legloc=0, ylim=None, labels=None, forced_single=False,
                    legend=False, xticks=[0.2, 0.4, 1.0, 2.0, 4.0, 10.0, 20.],
                    moleloc=None, eps=True):
    """

    Parameters
    ----------
    fname : str or list
        Path + filename of SMART transit file(s)
    """

    # Check fname
    if (type(fname) is str) or (type(fname) is unicode):
        # Single rad file provided
        multi = False
        # Make things lists
        fname = [fname]
        title = [title]
        labels = [labels]
        plot_kwargs = [plot_kwargs]
        # Construct figure and axes
        fig, ax = plt.subplots(figsize=(14,6))
    elif type(fname) is list:
        # Multiple rad files provided
        if labels is None: labels = [None for f in fname]
        if type(plot_kwargs) is not list: plot_kwargs = [plot_kwargs for f in fname]
        # Construct figure and axes
        if forced_single:
            multi = False
            fig, ax = plt.subplots(figsize=(14,6))
        else:
            multi = True
            fig, axs = plt.subplots(nrows=len(fname), figsize=(14,(6+1)*len(fname)))
    else:
        print "Invalid fname..."
        return

    # Loop over outputs
    for i in range(len(fname)):
        if multi:
            ax = axs[i]
        f = fname[i]
        label = labels[i]
        # Read in rad file(s)
        lam, zeff, tdepth = read_transit(f, radius=radius)
        mask = (lam > lammin) & (lam < lammax)
        # Set labels
        ax.set_xlabel(r"Wavelength [$\mu$m]")
        ax.set_ylabel(r"Atmosphere Radius [km]")
        ax.set_xlim([lammin, lammax])
        ax.semilogx()
        ax.text(0.5, 0.98, title[i],\
             verticalalignment='top', horizontalalignment='center',\
             transform=ax.transAxes,\
             color='black')
        if ylim is not None:
            ax.set_ylim([ylim[0], ylim[1]])
        # Actually plot it!
        ax.plot(lam[mask], zeff[mask], **plot_kwargs[i])

        if moleloc is not None:
            if type(moleloc) is list:
                mloc = moleloc[i]
            elif type(moleloc) is dict:
                mloc = moleloc
            else:
                print "Incompatible type for moleloc kwarg"
            if mloc is not None:
                for key, value in mloc.iteritems():
                    #print key, value
                    # get molecule color
                    mcol = molecules.color_from_molecule(key)
                    for im in range(len(value)):
                        # place label
                        ax.text(value[im][0], value[im][1], key, va='center', ha='center',
                             color=mcol, fontsize=16, zorder=10,
                             bbox=dict(boxstyle="square", fc="none", ec="none"))

        ax.set_xticks(xticks)
        ax.set_xticklabels(['{:g}'.format(tick) for tick in xticks])


        if multi or (i < 1):
            ax1 = ax.twinx()
        scale = 1e6
        ax1.plot(lam[mask], scale*tdepth[mask], visible=False, **plot_kwargs[i])
        ax1.set_ylim([scale*tdepth_from_zeff(ax.get_ylim()[0], radius[0], radius[1]),
                      scale*tdepth_from_zeff(ax.get_ylim()[1], radius[0], radius[1])])
        yticks = ax.get_yticks()
        ax1.set_ylabel(r"Transit Depth $(R_p/R_s)^2$ [ppm]", rotation=270, labelpad=30)
        #td_ticks = [tdepth_from_zeff(tick, radius[0], radius[0]) for tick in yticks]
        #ax1.set_yticks(td_ticks)
        #ax1.set_yticklabels(['{:g}'.format(age) for tick in ages.value]);

        if label is not None:
            # Get x and y data
            xdat = ax.lines[0].get_xdata()
            ydat = ax.lines[0].get_ydata()
            # Get x and y axis limits
            xmin, xmax = ax.axes.get_xlim()
            ymin, ymax = ax.axes.get_ylim()

            # degrade spec and put on even grid
            #dlam = (lammax - lammin) / 3000.
            #R = 10
            #xlow, dxlow = cg.noise_routines.construct_lam(lammin, lammax, dlam=dlam)
            #ylow = cg.downbin_spec(ydat, xdat, xlow, dlam=dxlow)
            #ax.plot(xlow, ylow, color="red")

            # Locate rolling max
            #df = pd.DataFrame(refl)
            df = pd.DataFrame(ydat)
            ymaxs = pd.rolling_max(df, int(len(ydat)/15.), min_periods=100)
            #ymaxs = ymaxs[mask].values.reshape(len(xdat))
            ymaxs = ymaxs.values.reshape(len(xdat))
            #ax.plot(xdat, ymaxs, color="blue")

            # Interpolate maxes to lower res
            dlam = (lammax - lammin) / 20.
            xlow, dxlow = cg.noise_routines.construct_lam(lammin, lammax, dlam=dlam)
            f2 = interp1d(xdat, ymaxs, kind='nearest', bounds_error=False, fill_value="extrapolate")
            ynew = f2(xlow)
            #ax.plot(xlow, ynew, color="green")

            # Interpolate back to high res grid
            f2 = interp1d(xlow, ynew, kind='nearest', bounds_error=False, fill_value="extrapolate")
            ycont = f2(xdat)

            # Extract relevant molecules
            molmask = np.array([True if molecules.molecules[j].keys()[0] in label else False
                                for j in range(len(molecules.molecules))])
            mols = molecules.molecules[molmask]
            # Set colors for each molecule
            colors = []
            for i in range(len(mols)):
                if i < 10:
                    ii = i
                else:
                    ii = i - 10
                colors.append("C%i" %ii)
            # Loop over relevant molecules
            for j in range(len(mols)):
                # Extract molecule formula (also dict key)
                key = mols[j].keys()[0]
                # Convert to latex text string
                text = molecules.tex_molecule(key)
                # Loop over bands in molecule
                for k in range(len(mols[j][key])):
                    # is bandcenter in plot window
                    if (mols[j][key][k] > lammin) & (mols[j][key][k] < lammax):
                        # Set x to be location of band
                        x = mols[j][key][k]
                        ix = np.argmin(np.fabs(xdat - x))
                        #y = np.random.uniform(yuse[ix], ymax)
                        y = ycont[ix] + np.random.uniform(0.01, 0.05)
                        #y = ydat[ix]
                        ax.text(x, y, text, va='center', ha='center',
                             color=colors[j], fontsize=12, zorder=100,
                             bbox=dict(boxstyle="square", fc="none", ec="none"))

        ax.set_xlim([lammin, lammax])

    # Label
    if legend:
        leg=ax.legend(loc=legloc, fontsize=16)
        leg.get_frame().set_alpha(0.0)

    # Save figure
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

    if eps:
        fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".eps"), bbox_inches='tight')

    return

def read_surf_rad(path, Nrow=14):
    """
    Parameters
    ----------
    path : str
        Path/name of file (*sur.rad)
    Nrow : int
        Number of elements in a row
    """


    # Convert each line to vector, compose array of vectors
    arrays = np.array([np.array(map(float, line.split())) for line in open(path)])

    # Flatten and reshape into rectangle grid
    arr = np.hstack(arrays).reshape((Nrow, -1), order='F')

    # Parse columns
    lam   = arr[0,:]
    wno   = arr[1,:]
    solar = arr[2,:]
    direct = arr[3,:]
    sky = arr[4,:]

    return lam, wno, solar, direct, sky

def plot_transmittance(fname, lammin=0.0, lammax=10.0, ylim=[0.0,1.0],
                       title="", plotdir="../../figures/", savetag="transmittance_test",
                       eps=True):
    """
    """

    # Read in file
    lam, wno, solar, direct, sky = read_surf_rad(fname)

    # Compute transmittance
    transmit = direct/solar

    # Make plot
    fig_params.set_spectrum_figsize()
    fig, ax = plt.subplots()
    ax.plot(lam, transmit, lw=2.0, color="black")
    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel(r"Transmittance")
    ax.set_title(title)
    ax.set_xlim(lammin, lammax)
    ax.set_ylim(ylim)

    # Save
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

    if eps:
        fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".eps"), bbox_inches='tight')
