"""
Figure 2

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    import numpy as np
    import os
    import sys
    sys.path.append("../../")
    from utils import fig_params
    import matplotlib.pyplot as plt
    from scipy.signal import savgol_filter

    # Params specific to this plot
    savetag = "fig02"
    plotdir="../../figures/"
    lammin  = 0.35
    lammax = 2.5

    files = [
        {
            "name" : "Conifer Forest",
            "file" : "conifers.txt",
            "header" : 26,
            "color" : "darkgreen",
            "ls" : "-"
        },
        {
            "name" : "Seawater",
            "file" : "seawater.txt",
            "header" : 26,
            "color" : "darkblue",
            "ls" : "-"
        },
        {
            "name" : "Ice/Snow",
            "file" : "snow.txt",
            "header" : 26,
            "color" : "lightblue",
            "ls" : "-"
        },
        {
            "name" : "Desert/Soil",
            "file" : "soil.txt",
            "header" : 26,
            "color" : "red",
            "ls" : "-"
        },
        {
            "name" : "Grassland",
            "file" : "grass.txt",
            "header" : 26,
            "color" : "lightgreen",
            "ls" : "-"
        },
        {
            "name" : "Goethite",
            "file" : "goethite.txt",
            "header" : 102,
            "color" : "darkred",
            "ls" : "-"
        },
        {
            "name" : "Basalt",
            "file" : "basalt_fresh.txt",
            "header" : 16,
            "color" : "grey",
            "ls" : "-"
        },
        {
            "name" : "Composite 1",
            "file" : "composite1.txt",
            "header" : 0,
            "color" : "black",
            "ls" : "dashed"
        },
        {
            "name" : "Composite 2",
            "file" : "composite2.txt",
            "header" : 0,
            "color" : "black",
            "ls" : "dotted"
        }
    ]

    #window1 = 11
    #window2 = 11
    #poly1 = 5
    #poly2 = 5

    fig, ax = plt.subplots(figsize=(10,6))
    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel(r"Geometric Albedo")
    ax.set_xlim([lammin, lammax])

    for i in range(len(files)):

        # Read in file
        path = os.path.join(os.path.dirname(__file__),"data_files/", files[i]["file"])
        data = np.genfromtxt(path, skip_header=files[i]["header"])
        wl, alb = data[:,0], data[:,1]

        mask = (wl > lammin) & (wl < lammax)

        # Plot
        ax.plot(wl[mask], alb[mask], color=files[i]["color"],
                                     ls=files[i]["ls"],
                                     label=files[i]["name"])



    # Label
    leg=ax.legend(loc=0, fontsize=16, ncol=2)
    leg.get_frame().set_alpha(0.0)

    # Save figure
    fig.savefig(os.path.join(os.path.dirname(__file__), plotdir, savetag+".pdf"), bbox_inches='tight')

    return
#########################################

if __name__ == "__main__":

    make_fig()
