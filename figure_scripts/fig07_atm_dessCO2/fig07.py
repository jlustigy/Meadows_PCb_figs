"""
Figure 7: Gao et al. atmospheric structure

Author: Jacob Lustig-Yaeger
"""

#########################################
def make_fig():

    import numpy as np
    from utils import fig_params
    import utils.atm as atm
    import os
    sys.path.append("../../")
    from utils import fig_params

    # Params specific to this plot
    savetag = "fig07"
    fname = "1barGao_final.atm"

    # Read in atm files
    atmpath = os.path.join(os.path.dirname(__file__),"model_outputs/", fname)
    data = np.genfromtxt(atmpath, skip_header=2)

    print data.shape
    # Feed to plotting function

    return
#########################################

if __name__ == "__main__":

    make_fig()
