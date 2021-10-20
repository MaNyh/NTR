import os
import numpy as np
import pyvista as pv
from tqdm import tqdm
import itertools

from NTR.database.case_dirstructure import casedirs
from NTR.utils.filehandling import yaml_dict_read
from NTR.utils.pyvista_utils import load_mesh, constructWallMesh, calc_dist_from_surface
from NTR.utils.geom_functions.distance import closest_node_index
from NTR.utils.mathfunctions import vecAbs


def calc_yplus(path_to_yaml_dict, verbose=True):
    """
    only suitable for structured meshes with orthogonal cells
    :param path_to_yaml_dict:
    :param verbose:
    :return:
    """
    settings = yaml_dict_read(path_to_yaml_dict)
    case_path = os.path.abspath(os.path.dirname(path_to_yaml_dict))

    use_velvar = settings["post_settings"]["yplus_calc"]["use_velfield"]
    use_rhofield = settings["post_settings"]["yplus_calc"]["use_rhofield"]

    volmesh_name = settings["post_settings"]["use_vtk_meshes"]["volmesh"]
    volmesh_path = os.path.join(case_path, casedirs["solution"], volmesh_name)
    wallnames = [os.path.join(case_path, casedirs["solution"], i) for i in
                            settings["post_settings"]["use_vtk_meshes"]["wallpatches"].values()]

    volmesh = load_mesh(volmesh_path)

    print("collecting walls...")
    walls = [load_mesh(w) for w in wallnames]
    walls_centers = [w.cell_centers() for w in walls]

    print("constructing near-wall-patches...")


    nearwalls = {}

    for idx in tqdm(range(volmesh.number_of_cells)):
        cell = volmesh.extract_cells(idx)
        edges = cell.extract_all_edges()
        for wallname, wall, wall_center in zip(wallnames, walls, walls_centers):

            found = False
            for edge_pt in edges.points:
                wallpt_idxs = closest_node_index(edge_pt, wall.points)
                wall_pt = wall.points[wallpt_idxs]
                if np.array_equal(edge_pt, wall_pt):
                    if wallname not in nearwalls.keys():
                        nearwalls[wallname] = cell
                    else:
                        nearwalls[wallname] = nearwalls[wallname].merge(cell)
                    found = True
                    break
            if found:
                break

    assert len(nearwalls.values())>0, "no connection from wall-meshes to the domain has been found"
    for nw_item, wall in zip(nearwalls.items(), walls):
        nw_name, nw_wall = nw_item
        nearwall_distancemap = calc_dist_from_surface(nw_wall.delaunay_2d(), wall)
        nearwall_dists = nw_wall.sample(nearwall_distancemap)
        nearwalls[nw_name].point_arrays["distances"] = nearwall_dists["distances"]
        nearwalls[nw_name] = nearwalls[nw_name].point_data_to_cell_data("distanes")

    for patchname, nearwall_mesh in nearwalls.items():
        print(patchname)
        #nearwall_mesh[use_velvar] = np.array(nearwall_vel)
        #nearwall_mesh["distances"] = np.array(center_ortho_walldist)
        nearwall_mesh["dudy"] = np.array([vecAbs(i) for i in nearwall_mesh[use_velvar]]) / np.array(
            nearwall_mesh["distances"])

        mu_0 = float(settings["simcase_settings"]["variables"]["DYNVISK"])

        if use_rhofield in nearwall_mesh.array_names:
            uTaus = getWalluTaus(mu_0, nearwall_mesh[use_rhofield], nearwall_mesh["dudy"])
        else:
            print("no rho found, assuming rho=1 (for incompressible flows set rho=1 and mu => nu)")
            uTaus = getWalluTaus(mu_0, np.ones(nearwall_mesh.number_of_cells), nearwall_mesh["dudy"])

        Deltay = nearwall_mesh["distances"] * uTaus / mu_0
        if verbose:
            nearwall_mesh["yplus"] = Deltay
            nearwall_mesh.set_active_scalars("yplus")
            nearwall_mesh.plot()
        print(wall)
        print("min : " + str(min(Deltay)))
        print("mean : " + str(np.mean(Deltay)))
        print("max : " + str(max(Deltay)))
    return 0 #Deltay


def getWalluTaus(mu_0, rhoW, gradUWall):
    tauW = gradUWall * mu_0 * rhoW
    u_tau = (tauW / rhoW) ** 0.5

    return u_tau
