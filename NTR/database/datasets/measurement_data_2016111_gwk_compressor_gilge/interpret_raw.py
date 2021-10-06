import numpy as np
import pyvista as pv
import os
from matplotlib import pyplot as plt
import pandas as pd

from NTR.preprocessing.create_geom import extract_geo_paras
from NTR.utils.mathfunctions import vecAbs
from NTR.utils.filehandling import write_pickle, read_pickle

if not os.path.isfile(("geometry.pkl")):
    print("no geometry.pkl found. will create one")

    pointfile = "profile_pointcloud.txt"
    profilepoints = np.loadtxt(pointfile) * 1e-3

    # =============================================================================
    # Bestimmung Profilparameter
    # =============================================================================
    sortedPoints, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, beta_meta_01, beta_meta_02, camber_angle = extract_geo_paras(
        profilepoints,
        alpha=0.01,
        verbose=False)

    x_mids = midsPoly.points[::, 0]
    y_mids = midsPoly.points[::, 1]
    x_ss = ssPoly.points[::, 0]
    y_ss = ssPoly.points[::, 1]
    x_ps = psPoly.points[::, 0]
    y_ps = psPoly.points[::, 1]

    stagger_angle = np.rad2deg(np.arctan((y_mids[-1] - y_mids[-0]) / (x_mids[-1] - x_mids[-0])))
    chordlength = vecAbs(midsPoly.points[0] - midsPoly.points[-1])

    geo_dict = {"points": sortedPoints,
                "sortedPoly": sortedPoints,
                "sidePolys": {"psPoly": psPoly, "ssPoly": ssPoly},
                "hk_vk_idx": {"ind_vk": ind_vk, "ind_hk": ind_hk},
                "midsPoly": midsPoly,
                "beta_metas": {"beta_meta_01": beta_meta_01, "beta_meta_02": beta_meta_02},
                "stagger_angle": stagger_angle,
                "camber_angle": camber_angle,
                }

    print("saving geometry.pkl")
    geo_filename = "geometry.pkl"
    write_pickle(os.path.join(os.path.dirname(__file__), geo_filename), geo_dict)

geo_dict = read_pickle("geometry.pkl")
print("geometry.pkl read")
camber = pv.Line((0, 0, 0), -(
        geo_dict["points"][geo_dict["hk_vk_idx"]["ind_vk"]] - geo_dict["points"][geo_dict["hk_vk_idx"]["ind_hk"]]))
xLine = pv.Line((-1, 0, 0), (1, 0, 0))
yLine = pv.Line((0, -1, 0), (0, 1, 0))

# this must be the test-data!
neu_pts = pv.PolyData(geo_dict["points"])
neu_pts.points -= neu_pts.points[geo_dict["hk_vk_idx"]["ind_vk"]]
neu_pts.rotate_z(-geo_dict["camber_angle"] + 90)

array_ = np.zeros(neu_pts.number_of_points)
camberlength = vecAbs(
    geo_dict["points"][geo_dict["hk_vk_idx"]["ind_vk"]] - geo_dict["points"][geo_dict["hk_vk_idx"]["ind_hk"]])

for idx, pt in enumerate(neu_pts.points):
    array_[idx] = pt[0] / camberlength

neu_pts["xc"] = array_

sortedPoints = geo_dict["points"]
"""
p = pv.Plotter()

p.add_mesh(camber)
p.add_mesh(xLine,color="black")
p.add_mesh(yLine)
p.add_mesh(sortedPoints)

p.add_mesh(neu_pts)
p.show()
"""


def calc_cp(px, pt1, p1):
    cp = (px - pt1) / (p1 - pt1)
    return cp


def read_gilgegwk():
    # global idx, ps_x, ps_cp, ss_x, ss_cp
    probe_camberposition_file = "raw_positions"
    probe_data_file = "raw_measuredata"
    position_arr = pd.read_csv(probe_camberposition_file, delimiter="\t")
    data_arr = pd.read_csv(probe_data_file, delimiter="\t")

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

    ps_x = []
    ps_p = []
    ps_cp = []

    ss_x = []
    ss_p = []
    ss_cp = []

    prel = pu
    pt1 = prel + popr
    p1 = prel + pspr

    nop_ps = 0
    nop_ss = 0
    for row in posvals:
        if row[1] == 1:
            ps_x.append(row[0] * 1e-3 / camberlength)
            nop_ps += 1
        elif row[2] == 1:
            ss_x.append(row[0] * 1e-3 / camberlength)
            nop_ss += 1
    for idx, p_mean_channel in enumerate(mean_channels.values()):
        if idx < nop_ps:
            ps_p.append(p_mean_channel + prel)
            ps_cp.append(calc_cp(ps_p[-1], pt1, p1))
        if idx >= nop_ps:
            ss_p.append(p_mean_channel + prel)
            ss_cp.append(calc_cp(ss_p[-1], pt1, p1))

    return ss_x, ss_cp, ps_x, ps_cp


ss_x, ss_cp, ps_x, ps_cp = read_gilgegwk()

fig, ax = plt.subplots()
# ToDo: Hier wird der plot "korrigiert". Das ist nicht schön. Wo stammt der Fehler her?
ax.plot(ps_x, -(np.array(ps_cp) - 1), "x", label="ps")
ax.plot(ss_x, -(np.array(ss_cp) - 1), "o", label="ss", )
ax.set(xlabel='x/c', ylabel='cp',
       title='profildruckverteilung gwk verdichter ')
ax.invert_yaxis()
ax.grid()
ax.legend()
plt.show()
