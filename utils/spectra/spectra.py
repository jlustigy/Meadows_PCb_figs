import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors

__all__ = ["read_rad", "plot_rad", "add_trans_plot"]

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
    ax.set_ylabel(r"Reflectivity ($I/F$)")

    ax.plot(lam, refl, **plot_kwargs)

    return

def plot_rad(fname, savetag="refl_test", lammin=0.2, lammax=20.0,
                    plot_kwargs={"color" : "black"}, title=None, plotdir="../../figures/",
                    legloc=None, ylim=None):
    """

    Parameters
    ----------
    fname : str or list
        Path + filename of SMART *.rad file
    """

    # Check fname
    if type(fname) is str:
        # Single rad file provided
        multi = False
        # Make things lists
        fname = [fname]
        title = [title]
        # Construct figure and axes
        fig, ax = plt.subplots(figsize=(14,6))
    elif type(fname) is list:
        # Multiple rad files provided
        multi = True
        # Construct figure and axes
        fig, axs = plt.subplots(nrows=len(fname), figsize=(14,(6+1)*len(fname)))
    else:
        print "Invalid fname..."
        return

    # Loop over outputs
    for i in range(len(fname)):
        if multi:
            ax = axs[i]
            f = fname[i]
        else:
            f = fname
        # Read in rad file(s)
        lam, wno, solar, toaf, rads = read_rad(f)
        mask = (lam > lammin) & (lam < lammax)
        # Calculate reflectivity
        refl = toaf / solar
        # Set labels
        ax.set_xlabel(r"Wavelength [$\mu$m]")
        ax.set_ylabel(r"Reflectivity ($I/F$)")
        ax.set_xlim([lammin, lammax])
        ax.text(0.5, 0.98, title[i],\
             verticalalignment='top', horizontalalignment='center',\
             transform=ax.transAxes,\
             color='black')
        if ylim is not None:
            ax.set_ylim([ylim[0], ylim[1]])
        # Actually plot it!
        ax.plot(lam[mask], refl[mask], **plot_kwargs)



    # Save figure
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

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
