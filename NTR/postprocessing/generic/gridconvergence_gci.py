import pyvista as pv
import numpy as np


def pEq(x, epsilon32, epsilon21, r21, r32, s):
    return ((np.log(epsilon32 / epsilon21) + np.log(((r21 ** x) - s) / ((r32 ** x) - s))) / np.log(r21) - x)


def getGCI(numcell_mesh1, numcell_mesh2, numcell_mesh3,fc1,fc2,fc3,grid_dimension, Fs=1.25):

    if ((fc1 > fc2) and (fc2 < fc3)):
        print('Error: fci not monoton! ')
        print('Check your input, repeat your grid convergence study, or contact your CFD-expert')
        return -1

    if ((fc1 < fc2) and (fc2 > fc3)):
        print('Error: fci not monoton!')
        print('Check your input, repeat your grid convergence study, or contact your CFD-expert')
        return -1

    if ((fc1 == fc2) or (fc2 == fc3)):
        print('Error: fci not monoton! ')
        print('Check your input, repeat your grid convergence study, or contact your CFD-expert')
        return -1

    if ((numcell_mesh1 <= numcell_mesh2) or (numcell_mesh2 <= numcell_mesh3)):
        print('Error: Number of Grid Nodes illogical! Check your input!')
        return -1

    print("input data valid")

    r21 = (numcell_mesh1 / numcell_mesh2) ** (1 / grid_dimension)
    r32 = (numcell_mesh2 / numcell_mesh3) ** (1 / grid_dimension)

    epsilon32 = (fc3 - fc2)
    epsilon21 = (fc2 - fc1)

    s = 1 * np.sign(epsilon32 / epsilon21)

    x0 = np.log(epsilon32 / epsilon21) / np.log(r21)
    p = pEq(x0, epsilon32, epsilon21, r21, r32, s)

    GCI_1 = Fs * abs((epsilon21 / fc1)) * 1 / ((r21 ** p) - 1)
    GCI_2 = Fs * abs((epsilon32 / fc2)) * 1 / ((r32 ** p) - 1)
    GCI_3 = GCI_2 * r32 ** p

    f_extra = fc1 + (fc1 - fc2) / ((r21 ** p) - 1)

    EERE_1 = abs((f_extra - fc1) / f_extra)
    EERE_2 = abs((f_extra - fc2) / f_extra)
    EERE_3 = abs((f_extra - fc3) / f_extra)

    A_Flag = GCI_2 / (GCI_1 * (r21 ** p))

    if ((A_Flag > 1.15) | (A_Flag < 0.85)):
        print('Error: Your fc-values are not in the asymptotic range!')
        print('Refine your grid and repeat grid convergence study, or contact your CFD-expert')

    GCI_1_p1 = Fs * abs((epsilon21 / fc1)) * 1 / ((r21) - 1)
    GCI_2_p1 = Fs * abs((epsilon32 / fc2)) * 1 / ((r32) - 1)
    GCI_3_p1 = GCI_2_p1 * r32

    f_extra_p1 = fc1 + (fc1 - fc2) / ((r21) - 1)

    EERE_1_p1 = abs((f_extra_p1 - fc1) / f_extra_p1)
    EERE_2_p1 = abs((f_extra_p1 - fc2) / f_extra_p1)
    EERE_3_p1 = abs((f_extra_p1 - fc3) / f_extra_p1)

    return [GCI_1, GCI_2, GCI_3], [EERE_1, EERE_2, EERE_3], [GCI_1_p1, GCI_2_p1, GCI_3_p1], \
           [EERE_1_p1, EERE_2_p1, EERE_3_p1]

mesh3_ncell = 400000
mesh2_ncell = 1200000
mesh1_ncell = 3600000

fc3 = 23.263
fc2 = 23.165
fc1 = 23.151

grid_dimension = 3

a = getGCI(mesh1_ncell, mesh2_ncell, mesh3_ncell,fc1,fc2,fc3,grid_dimension, Fs=1.25)
