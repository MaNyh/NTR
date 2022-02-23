import pyvista as pv
import numpy as np
from tqdm import tqdm
import os

from NTR.utils.mesh_handling.pyvista_utils import load_mesh
from NTR.utils.filehandling import yaml_dict_read
from NTR.database.case_dirstructure import casedirs
from NTR.utils.mathfunctions import vecAbs, lineseg_dist


def vol_to_line(vtkmesh, ave_direction, verbose=False):
    """
    this function is assuming a structured grid without curved gridlines
    it extracts layers and averages them. currently the face-normals have to be alligned with the global coordinate system

    :param settings:
    :param verbose:
    :return:
    """
    mesh = vtkmesh
    array_names_raw = mesh.array_names
    array_names = []
    for key in array_names_raw:
        if key not in array_names:
            array_names.append(key)

    dirs = {"x": 0, "y": 2, "z": 4}
    interpol_dir = dirs[ave_direction]

    rest = mesh.copy()

    pts = []
    meanvals = {}
    for array_name in array_names:
        meanvals[array_name] = []

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

        centers_poi = centers.points[::, int(interpol_dir - interpol_dir / 2)]
        points_bnd = np.ones(len(centers.points)) * bnd

        ids = np.where(np.equal(centers_poi, points_bnd))[0]
        ids_negative = np.where(np.not_equal(centers_poi, points_bnd))[0]

        assert mesh.number_of_cells == (len(ids) + len(ids_negative) + mesh.number_of_cells - rest.number_of_cells), \
            "somethings wrong"

        layer = rest.extract_cells(ids)

        if len(ids_negative) > 0:
            rest = rest.extract_cells(np.array([i for i in range(len(centers.points)) if not np.isin(i, ids)]))
        else:
            rest = pv.UniformGrid()

        for array_name in array_names:
            mean = layer[array_name].mean(axis=0)
            meanvals[array_name].append(mean)

        pts.append(bnd)
        pbar.update(len(ids))
    pbar.close()
    pos = np.array(pts)
    vals = {}
    for array_name in array_names:
        vals[array_name] = np.array(meanvals[array_name])

    return pos, vals



def vol_to_plane(mesh, ave_direction, verbose=False):
    """
    this function is assuming an extruded mesh in the average-direction.
    it extracts cells along the average direction and it is averaging all arrays of this cell along the extracted set of cells

    :param mesh: a pyvista-mesh
    :param ave_direction: char (x,y,z)
    :param verbose: show progress in plots
    :return: merged averaged cells into pv.PolyData-format
    """

    array_names = mesh.array_names

    dirs = {"x": 0, "y": 2, "z": 4}
    interpol_dir = dirs[ave_direction]
    rest = mesh.copy()


    pbar = tqdm(total=mesh.number_of_cells)
    slices = []
    while (rest.number_of_cells > 0):
        if verbose:
            p = pv.Plotter()
            p.add_mesh(mesh, opacity=0.5)
            p.add_mesh(rest)
            p.show()

        bounds = list(rest.bounds)
        bounds_centers = list(rest.cell_centers().bounds)

        bounds[interpol_dir] = bounds_centers[interpol_dir+1]


        cells_on_line_ids = rest.find_cells_within_bounds(bounds)
        ids_negative = [i for i in range(rest.number_of_cells) if i not in cells_on_line_ids]

        assert mesh.number_of_cells == (len(cells_on_line_ids) + len(ids_negative) + mesh.number_of_cells - rest.number_of_cells), \
            "somethings wrong"

        slice_ave = rest.extract_cells(cells_on_line_ids)

        slices.append(slice_ave)

        if len(ids_negative) > 0:
            rest = rest.extract_cells(
                np.array([i for i in range(rest.number_of_cells) if not np.isin(i, cells_on_line_ids)]))
        else:
            rest = pv.UniformGrid()

        pbar.update(len(cells_on_line_ids))
    slice_ave = slices[0].copy()
    slice_ave.clear_data()
    for s in slices:
        for arname in array_names:
            if arname not in slice_ave.array_names:
                slice_ave.point_data[arname] = s.point_data[arname]
            else:
                slice_ave.point_data[arname] += s.point_data[arname]
    for arname in array_names:
        slice_ave.point_data[arname]/=len(slices)
    slice_ave.point_data_to_cell_data()
    pbar.close()


    return slice_ave


def vol_to_plane_fromsettings(settings_yml_path):
    settings = yaml_dict_read(settings_yml_path)
    casepath = os.path.abspath(os.path.dirname(settings_yml_path))
    meshpath = os.path.join(casepath, casedirs["solution"], settings["post_settings"]["use_vtk_meshes"]["volmesh"])
    line_direction = settings["post_settings"]["average_volumeonplane"]["line_dir"]
    cellcentered = settings["post_settings"]["average_volumeonplane"]["cellcentered"]
    mesh = load_mesh(meshpath)

    ave_slice = vol_to_plane(mesh, line_direction, cell_centered=cellcentered)
    return ave_slice


def vol_to_line_fromsettings(settings_yml_path):
    settings = yaml_dict_read(settings_yml_path)
    casepath = os.path.abspath(os.path.dirname(settings_yml_path))
    meshpath = os.path.join(casepath, casedirs["solution"], settings["post_settings"]["use_vtk_meshes"]["volmesh"])
    line_direction = settings["post_settings"]["average_volumeonline"]["line_dir"]

    mesh = load_mesh(meshpath)
    points, data = vol_to_line(mesh, line_direction)
    return points, data
