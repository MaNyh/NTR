import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

from NTR.utils.pyvista_utils import load_mesh
from NTR.utils.filehandling import yaml_dict_read
from NTR.database.case_dirstructure import casedirs

def vol_to_line_fromsettings(settings_yml_path):
    settings = yaml_dict_read(settings_yml_path)
    casepath = os.path.abspath(os.path.dirname(settings_yml_path))
    meshpath = os.path.join(casepath,casedirs["solution"],settings["post_settings"]["use_vtk_meshes"]["volmesh"])
    line_direction = settings["post_settings"]["average_volumeonline"]["line_dir"]

    mesh = load_mesh(meshpath)
    ans = vol_to_line(mesh,line_direction)
    print(ans)

def vol_to_line(vtkmesh, ave_direction, verbose=False):
    """
    this function is assuming a structured grid without curved gridlines
    it extracts layers and averages them. currently the face-normals have to be alligned with the global coordinate system

    :param settings:
    :param verbose:
    :return:
    """
    mesh = vtkmesh
    array_names = mesh.array_names

    dirs = {"x": 0, "y": 2, "z": 4}
    interpol_dir = dirs[ave_direction]

    rest = mesh.copy()
    pts = []
    meanvals = {}
    pbar = tqdm(total=mesh.number_of_cells)

    while (rest.number_of_cells > 0):
        if verbose:
            p = pv.Plotter()
            p.add_mesh(mesh, opacity=0.5)
            p.add_mesh(rest)
            p.show()
        centers = rest.cell_centers()
        bounds = centers.bounds
        bnd = bounds[interpol_dir]

        ids = np.where(
            np.equal(centers.points[::, int(interpol_dir - interpol_dir / 2)], np.ones(len(centers.points)) * bnd))[0]

        ids_negative = np.where(
            np.not_equal(centers.points[::, int(interpol_dir - interpol_dir / 2)], np.ones(len(centers.points)) * bnd))[0]

        if len(ids_negative) > 0:
            rest = rest.extract_cells(np.array([i for i in range(len(centers.points)) if not np.isin(i, ids)]))
        else:
            rest = pv.UniformGrid()
        layer = mesh.extract_cells(ids)

        for array_name in array_names:
            meanvals[array_name] = []

        for array_name in array_names:
            mean = np.average(layer[array_name], axis=0)
            meanvals[array_name].append(mean)

        pts.append(bnd)
        pbar.update(len(ids))
    pbar.close()
    pos = np.array(pts)
    vals = {}
    for array_name in array_names:
        vals[array_name] = nnp.array(meanvals[array_name])

    return pos, vals
