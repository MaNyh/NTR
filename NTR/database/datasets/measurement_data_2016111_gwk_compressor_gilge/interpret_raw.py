import numpy as np
import pyvista as pv
import os
from matplotlib import pyplot as plt
import pandas as pd

from NTR.utils.geom_functions.profileparas import calcMidPassageStreamLine
from NTR.preprocessing.create_geom import extract_geo_paras
from NTR.utils.mathfunctions import vecAbs
from NTR.utils.filehandling import write_pickle, read_pickle

if not os.path.isfile(("geometry.pkl")):
    print("no geometry.pkl found. will create one")

    pointfile = "profile_pointcloud.txt"
    profilepoints = np.loadtxt(pointfile)*1e-3

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
                "camber_angle":camber_angle,
                }

    print("saving geometry.pkl")
    geo_filename = "geometry.pkl"
    write_pickle(os.path.join(os.path.dirname(__file__), geo_filename), geo_dict)

geo_dict = read_pickle("geometry.pkl")
print("geometry.pkl read")
camber = pv.Line((0,0,0),-(geo_dict["points"][geo_dict["hk_vk_idx"]["ind_vk"]]-geo_dict["points"][geo_dict["hk_vk_idx"]["ind_hk"]]))
xLine = pv.Line((-1,0,0),(1,0,0))
yLine = pv.Line((0,-1,0),(0,1,0))

#this must be the test-data!
neu_pts = pv.PolyData(geo_dict["points"])
neu_pts.points-=neu_pts.points[geo_dict["hk_vk_idx"]["ind_vk"]]
neu_pts.rotate_z(-geo_dict["camber_angle"] + 90)

array_ = np.zeros(neu_pts.number_of_points)
camberlength = vecAbs(geo_dict["points"][geo_dict["hk_vk_idx"]["ind_vk"]]-geo_dict["points"][geo_dict["hk_vk_idx"]["ind_hk"]])

for idx,pt in enumerate(neu_pts.points):
    array_[idx]=pt[0]/camberlength

neu_pts["xc"]=array_

sortedPoints=geo_dict["points"]

p = pv.Plotter()

p.add_mesh(camber)
p.add_mesh(xLine,color="black")
p.add_mesh(yLine)
p.add_mesh(sortedPoints)

p.add_mesh(neu_pts)
p.show()

probe_camberposition_file = "raw_positions"
probe_data_file = "raw_measuredata"

position_arr = pd.read_csv(probe_camberposition_file,delimiter="\t")
data_arr = pd.read_csv(probe_data_file,delimiter="\t")
