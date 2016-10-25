# Executing this script will generate all figures for the manuscript Meadows et al. 2016 (submitted).

import os
import numpy as np
import imp

##########################################
def filter_hidden(listdir):
    # Filter out hidden files: .*
    igood = []
    for i in range(len(listdir)):
        if not listdir[i].startswith("."):
            igood.append(i)
    listdir = np.array(listdir)
    return list(listdir[np.array(igood, dtype=int)])

def filter_files(listdir):
    # Filter out files with "."
    igood = []
    for i in range(len(listdir)):
        if "." not in listdir[i]:
            igood.append(i)
    listdir = np.array(listdir)
    return list(listdir[np.array(igood, dtype=int)])

def generate_all_figures(topdir="figure_scripts/"):

        # Get names of all figure script directories
        listdir = filter_hidden(os.listdir(topdir))

        listdir = filter_files(listdir)

        # Grab all python files and store them in a list
        pyfiles = []
        # Loop over all fig dirs
        for i in range(len(listdir)):
            listeach = filter_hidden(os.listdir(os.path.join(topdir,listdir[i])))
            # Loop over all files within
            for j in range(len(listeach)):
                # If it's a python file add full path to the list
                if listeach[j].endswith(".py"):
                    pyfiles.append(os.path.join(topdir,listdir[i],listeach[j]))

        for path in pyfiles:
            foo = imp.load_source('module.name', path)
            foo.make_fig()

##########################################
if __name__ == "__main__":

    generate_all_figures()
