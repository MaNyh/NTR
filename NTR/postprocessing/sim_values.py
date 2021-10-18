import os
import numpy as np

from NTR.database.case_dirstructure import casedirs
from NTR.utils.filehandling import yaml_dict_read
from NTR.utils.pyvista_utils import load_mesh, constructWallMesh
from NTR.utils.geom_functions.distance import closest_node_index
from NTR.utils.mathfunctions import vecAbs, vecDir


def calc_yplus(path_to_yaml_dict,verbose=True):
    settings = yaml_dict_read(path_to_yaml_dict)
    case_path = os.path.abspath(os.path.dirname(path_to_yaml_dict))

    use_velvar = settings["post_settings"]["yplus_calc"]["use_velfield"]
    use_rhofield = settings["post_settings"]["yplus_calc"]["use_rhofield"]

    volmesh_name = settings["post_settings"]["use_vtk_meshes"]["volmesh"]
    volmesh_path = os.path.join(case_path, casedirs["solution"], volmesh_name)
    wallpatches_namelist = [os.path.join(case_path, casedirs["solution"], i) for i in
                            settings["post_settings"]["use_vtk_meshes"]["wallpatches"].values()]

    volmesh = load_mesh(volmesh_path)
    volmesh_centers = volmesh.cell_centers()
    wall = constructWallMesh(wallpatches_namelist)
    wall_centers = wall.cell_centers().points

    print("extracting near-wall cells...")
    nearwall_ids = []
    nearwall_vel = []
    center_ortho_walldist = []

    for idx, cellcenter in enumerate(wall_centers):
        nearwall_idx = closest_node_index(cellcenter, volmesh_centers.points)
        nearwall_ids.append(nearwall_idx)
        nearwall_vel.append(volmesh_centers[use_velvar][nearwall_idx])
        center_ortho_walldist.append(vecAbs(cellcenter - volmesh_centers.points[nearwall_idx]))

    nearwall_mesh = volmesh.extract_cells(nearwall_ids)
    nearwall_mesh[use_velvar] = np.array(nearwall_vel)
    nearwall_mesh["dist"] = np.array(center_ortho_walldist)
    nearwall_mesh["dudy"] = np.array([vecAbs(i) for i in nearwall_mesh[use_velvar]])

    mu_0 = float(settings["simcase_settings"]["variables"]["DYNVISK"])

    uTaus = getWalluTaus(mu_0, nearwall_mesh[use_rhofield],nearwall_mesh["dudy"])
    Deltay = nearwall_mesh["dist"] * uTaus / mu_0
    if verbose:
        nearwall_mesh["yplus"]=Deltay
        nearwall_mesh.set_active_scalars("yplus")
        nearwall_mesh.plot()
    return Deltay


def getWalluTaus( mu_0, rhoW, gradUWall):

    tauW = gradUWall * mu_0 * rhoW
    u_tau = (tauW / rhoW) ** 0.5

    return u_tau
