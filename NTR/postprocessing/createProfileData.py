import matplotlib.pyplot as plt
import os
import pyvista as pv
import numpy as np
from NTR.utils.geom_functions import GetProfileValuesMidspan


def createProfileData(path_to_mesh):
    [value_names, [values_ss, values_ps]] = GetProfileValuesMidspan(path_to_mesh)

    x_ss = list(values_ss[value_names.index('X')])[::-1]
    y_ss = list(values_ss[value_names.index('Y')])[::-1]
    p_ss = list(values_ss[value_names.index('p')])[::-1]

    x_ps = list(values_ps[value_names.index('X')])
    y_ps = list(values_ps[value_names.index('Y')])
    p_ps = list(values_ps[value_names.index('p')])

    plt.figure(figsize=(8, 8))
    plt.plot(x_ss, y_ss, '-r', lw=1)
    plt.plot(x_ps, y_ps, '-b', lw=1)

    plt.axis('equal')
    output_path = os.path.dirname(os.path.abspath(path_to_mesh))
    plt.grid()
    plt.savefig(os.path.join(output_path, 'kontrollplot_profil.pdf'))
