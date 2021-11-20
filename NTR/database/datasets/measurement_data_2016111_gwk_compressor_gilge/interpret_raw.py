import os.path

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os

from NTR.postprocessing.profile_loading import calc_inflow_cp

def read_gilgegwk(verbose=True):

    probe_camberposition_file = "raw_positions"
    probe_data_file = "raw_measuredata"
    basedir = os.path.dirname(__file__)
    position_fname = os.path.join(basedir,probe_camberposition_file)
    position_arr = pd.read_csv(position_fname, delimiter="\t")
    data_fname = os.path.join(basedir,probe_data_file)
    data_arr = pd.read_csv(data_fname, delimiter="\t")
    #todo: where get this data from? hardcoding no good!
    camberlength = 0.13
    posvals = [[float(f) for f in i[0].split()] for i in position_arr.values[:, :-2]]
    head = data_arr.keys()
    mean_channels = {}

    for idx, channel in enumerate(head[2:13]):
        mean_channels[channel + "_mean"] = np.mean([float(i) for i in data_arr.loc[1:, channel].values])
        posvals[idx].append(mean_channels[channel + "_mean"])

    # mach = np.mean([float(i) for i in data_arr.loc[1:, "Machzahl"].values])
    pu = np.mean([float(i) for i in data_arr.loc[1:, "Pu"].values])
    popr = np.mean([float(i) for i in data_arr.loc[1:, "PoPr"].values])
    pspr = np.mean([float(i) for i in data_arr.loc[1:, "PsPr"].values])

    ps_xc = []
    ps_p = []
    ps_cp = []

    ss_xc = []
    ss_p = []
    ss_cp = []

    prel = pu
    pt1 = prel + popr
    p1 = prel + pspr

    nop_ps = 0
    nop_ss = 0
    for row in posvals:
        if row[1] == 1:
            ps_xc.append(row[0] * 1e-3 / camberlength)
            nop_ps += 1
        elif row[2] == 1:
            ss_xc.append(row[0] * 1e-3 / camberlength)
            nop_ss += 1
    for idx, p_mean_channel in enumerate(mean_channels.values()):
        if idx < nop_ps:
            ps_p.append(p_mean_channel + prel)
            ps_cp.append(calc_inflow_cp(ps_p[-1], pt1, p1))
        if idx >= nop_ps:
            ss_p.append(p_mean_channel + prel)
            ss_cp.append(calc_inflow_cp(ss_p[-1], pt1, p1))

    if verbose:
        fig, ax = plt.subplots()
        # ToDo: Hier wird der plot "korrigiert". Das ist nicht schön. Wo stammt der Fehler her?
        ax.plot(ps_xc, -(np.array(ps_cp) - 1), "x", label="ps")
        ax.plot(ss_xc, -(np.array(ss_cp) - 1), "o", label="ss", )
        ax.set(xlabel='x/c', ylabel='cp',
               title='profildruckverteilung gwk verdichter ')
        ax.invert_yaxis()
        ax.grid()
        ax.legend()
        plt.show()

    return ss_xc, ss_cp, ps_xc, ps_cp


