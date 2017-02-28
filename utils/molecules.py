
import numpy as np

def colors_from_molecs(molecules):
    return

def tex_molecule(formula):
    """
    Create LaTeX compatible string for molecular formula from simple string
    """
    new = r''
    for char in formula:
        if char.isdigit():
            tmp = '$_{'+char.upper()+'}$'
        else:
            tmp = char
        new = new+tmp
    return new

molecules = np.array([
    {"O2" : [0.2, 0.69, 0.76, 1.27]},
    {"O3" : [2.5, 0.55, 9.6]},
    {"O4" : [0.45, 0.48, 0.53, 0.57, 0.63, 1.06, 1.27]},
    {"CH4" : [0.79, 0.89, 1.0, 1.1, 1.4, 1.7, 2.31, 3.3, 7.7]},
    {"CO2" : [1.05, 1.21, 1.6, 2.01, 4.2, 9.4, 10.4, 15.0]},
    {"CO" : [1.6, 2.35]},
    {"N4" : [4.1]},
    {"N2O" : [2.11, 2.25, 2.6, 2.67, 2.97, 3.6, 3.9, 4.3, 4.5, 7.9, 17.0]},
    {"H2" : [0.65, 0.825]},
    {"H2O" : [0.65, 0.72, 0.82, 0.94, 1.12, 1.4, 1.85, 2.5, 6.3]}
])
