import os
import numpy as np
from scipy import integrate
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors

__all__ = ["plot_single_atm", "add_atm_plot", "plot_double_atm", "plot_aerosol"]

def add_atm_plot(P, T, gas_profiles, molec_names, ax0=None, legend=False,
                 title=None, xlim=None, ylim=None, tlim=None, bbox_to_anchor=(1., 1.03),
                 seed=0, legloc=None):
    """
    """

    # Create new figure and axis or use the one provided
    if ax0 is None:
        fig, ax = plt.subplots(figsize=(7,8))
    else:
        ax = ax0

    axT = ax.twiny()
    ax.invert_yaxis()
    ax.loglog()
    ax.set_xlabel("Volume Mixing Ratio")
    ax.set_ylabel("Pressure [Pa]")
    axT.set_xlabel("Temperature [K]")
    axT.set_axisbelow(True)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if tlim is not None:
        axT.set_xlim(tlim)

    ax.plot(0,0,color="black", label="Temp.", ls="--")

    colors = []
    for i in range(len(molec_names)):
        if i < 10:
            ii = i
        else:
            ii = i - 10
        colors.append("C%i" %ii)

    # Plot gases
    for i in range(len(molec_names)):
        ax.plot(gas_profiles[:,i], P, label=molec_names[i], color=colors[i])

    # Plot temperature structure
    axT.plot(T, P, color="black", label="Temperature", ls="--")

    # Create legend
    if legend:
        if legend == "custom":
            # Set random seed
            np.random.seed(seed)
            # Get axis limits
            xmin, xmax = ax.axes.get_xlim()
            tmin, tmax = axT.axes.get_xlim()
            ymin, ymax = ax.axes.get_ylim()
            # Create masks including only visible data
            tmask = (T > (tmin+tmin/20)) & (P < (tmin+tmin/20))
            ymask = (P > (ymax*2)) & (P < (ymin/2))
            mask = (tmask) & (ymask)
            # Does temperature line show up in the plot window?
            if np.sum(mask) > 0:
                # Select y-index for text label
                if legloc is None:
                    # Grab a valid y-index at random
                    j = np.random.choice(np.arange(len(mask))[mask])
                elif legloc[0] is None:
                    # Grab a valid y-index at random
                    j = np.random.choice(np.arange(len(mask))[mask])
                else:
                    # Grab y-index nearest user specified y-value
                    j = np.argmin(np.fabs(P - legloc[0]))
                # Set x, y of text label
                y, x = P[j], T[j]
                # Place text
                axT.text(x, y, "Temp", va='center', ha='center',
                    color="black", fontsize=12, zorder=100,
                    bbox=dict(boxstyle="square", fc="w", ec="none"))
            # Loop over molecules
            for i in range(len(molec_names)):
                xmask = (gas_profiles[:,i] > (xmin*2)) & (gas_profiles[:,i] < (xmax/2))
                mask = (xmask) & (ymask)
                # Does the line show up in the plot window?
                if np.sum(mask) > 0:
                    # Select y-index for text label
                    if legloc is None:
                        # Grab a valid y-index at random
                        j = np.random.choice(np.arange(len(mask))[mask])
                    elif legloc[i+1] is None:
                        # Grab a valid y-index at random
                        j = np.random.choice(np.arange(len(mask))[mask])
                    else:
                        # Grab y-index nearest user specified y-value
                        j = np.argmin(np.fabs(P - legloc[i+1]))
                    # Set x, y of text label
                    y, x = P[j], gas_profiles[j,i]
                    # Place text
                    ax.text(x, y, molec_names[i], va='center', ha='center',
                         color=colors[i], fontsize=12, zorder=100,
                         bbox=dict(boxstyle="square", fc="w", ec="none"))
                    """
                    ns = np.arange(len(mask))*mask
                    nmin, nmax = np.min(np.nonzero(ns)), ns.max()
                    #jxmin = np.argmax(np.fabs(gas_profiles[mask,i] - xmin))
                    #jxmax = np.argmax(np.fabs(gas_profiles[mask,i] - xmax))
                    #xmindiff = np.power(np.log10(gas_profiles[mask,i]) - np.log10(xmin), 2)
                    #xmaxdiff = np.power(np.log10(gas_profiles[mask,i]) - np.log10(xmax), 2)
                    #ymindiff = np.power(np.log10(P[mask]) - np.log10(ymin), 2)
                    #ymaxdiff = np.power(np.log10(P[mask]) - np.log10(ymax), 2)
                    """
                    """
                    if np.sum(xmindiff) > np.sum(xmaxdiff):
                        xdiff = xmaxdiff
                    else:
                        xdiff = xmindiff
                    if np.sum(ymindiff) > np.sum(ymaxdiff):
                        ydiff = ymaxdiff
                    else:
                        ydiff = ymindiff
                    """
                    """
                    #vmin, vmax = np.min(gas_profiles[mask,i]), np.max(gas_profiles[mask,i])
                    #j = np.argmin(gas_profiles[mask,i] - np.power(10, 0.5*(np.log10(vmax) + np.log10(vmin))))
                    #lspace =  (np.log10(np.hstack([gas_profiles[mask, :i], gas_profiles[mask,i+1:]])) - np.log10(gas_profiles[mask,i, np.newaxis]))
                    #jj = np.argmin(np.sum(lspace, axis=0)) # Theoretically the closest line
                    #j = np.argmin(xmindiff+xmaxdiff+lspace[:,jj])
                    #j = np.argmax(lspace[:,jj])
                    #j = np.argmax(xdiff)
                    j = np.random.randint(nmin, nmax)
                    j = np.random.choice(np.arange(len(mask))[mask])
                    y, x = P[j], gas_profiles[j,i]
                    #y, x = P[mask][j], gas_profiles[mask,i][j]
                    ax.text(x, y, molec_names[i], va='center', ha='center',
                         color=colors[i], fontsize=12, zorder=100,
                         bbox=dict(boxstyle="square", fc="w", ec="none"))
                    """
        else:
            # Regular legend
            leg=ax.legend(fontsize=12, bbox_to_anchor=bbox_to_anchor)
            leg.get_frame().set_alpha(0.0)
    else:
        # No legend
        pass

    # Add title
    if title is not None:
        title = r"\underline{%s}" %title
        ax.set_title(title, y=1.1)

    if ax0 is None:
        #plt.show()
        fig.savefig("test.pdf", bbox_inches="tight")
        return
    else:
        return


def plot_single_atm(P, T, gas_profiles, molec_names, savetag="atm_test",
                    title=None, plotdir="../../figures/", bbox_to_anchor=(1., 1.03),
                    xlim=None, ylim=None, tlim=None, seed=0, legend="custom",
                    legloc=None):
    """
    """

    # Create figure
    fig, ax = plt.subplots(figsize=(7,8))

    # Add atm plot
    add_atm_plot(P, T, gas_profiles, molec_names, ax0=ax, legend=legend,
                     title=title, bbox_to_anchor=bbox_to_anchor, xlim=xlim,
                     ylim=ylim, tlim=tlim, seed=seed, legloc=legloc)

    # Save figure
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

def plot_double_atm(P, T, gas_profiles, molec_names, savetag="atm_test",
                    title=None, plotdir="../../figures/", bbox_to_anchor=(1., 1.03),
                    xlim=None, ylim=None, tlim=None, legend="custom", seed=0,
                    legloc=None):
    """
    """

    # Create figure
    fig, ax = plt.subplots(1,2, figsize=(15,8))

    # Add atm plot
    add_atm_plot(P[0], T[0], gas_profiles[0], molec_names[0], ax0=ax[0], legend=legend,
                     title=title[0], bbox_to_anchor=bbox_to_anchor, xlim=xlim, ylim=ylim,
                     tlim=tlim, seed=seed, legloc=legloc)
    add_atm_plot(P[1], T[1], gas_profiles[1], molec_names[1], ax0=ax[1], legend=legend,
                     title=title[1], xlim=xlim, ylim=ylim, tlim=tlim, seed=seed,
                     legloc=legloc)

    # Save figure
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

    return

def plot_aerosol(altitude, density, radius, savetag="aero_test",
                    plotdir="../../figures/"):
    """
    """

    # Create figure
    fig, ax = plt.subplots(figsize=(7,8))
    ax2 = ax.twiny()
    ax.semilogx()
    ax2.semilogx()
    ax.set_xlabel(r"Aerosol Number Density [cm$^{-3}$]")
    ax2.set_xlabel(r"Particle Radius [$\mu$m]", labelpad=10)
    ax.set_ylabel(r"Altitude [km]")
    ax.set_ylim([0,100])
    ax.set_xlim([1e-1, 1e5])
    ax2.set_xlim([1e-3, 1e0])

    # Add atm plot
    #import pdb; pdb.set_trace()
    ax2.plot(radius, altitude, label="Radius", color="black", ls="--")
    ax.plot(1, 1, label="Radius", color="black", ls="--")
    ax.plot(density, altitude, label="Density", color="black")

    # Label
    leg=ax.legend(loc=6, fontsize=16)
    leg.get_frame().set_alpha(0.0)

    # Save figure
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

    return
