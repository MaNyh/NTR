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


def vol_to_plane(volmesh, ave_direction, cell_centered=False, verbose=False):
    volume = volmesh
    if cell_centered:
        cell_centers = volume.cell_centers()
        mesh = cell_centers
    else:
        mesh = volume
    dirs = {"x": 0, "y": 2, "z": 4}
    interpol_dir = dirs[ave_direction]
    boundshigh = mesh.bounds[interpol_dir]
    boundslow = mesh.bounds[interpol_dir + 1]

    helper_one = (interpol_dir + 2) % 6
    helper_one_low = mesh.bounds[helper_one + 1]

    helper_two = (interpol_dir + 4) % 6
    helper_two_low = mesh.bounds[helper_two + 1]

    end = [None, None, None]
    end[int(interpol_dir / 2)] = boundslow
    end[int(helper_one / 2)] = helper_one_low
    end[int(helper_two / 2)] = helper_two_low
    end = np.array(end)

    base = [None, None, None]
    base[int(interpol_dir / 2)] = boundshigh
    base[int(helper_one / 2)] = helper_one_low
    base[int(helper_two / 2)] = helper_two_low
    base = np.array(base)

    pts = []
    tolerance = vecAbs(base - end) / 1000
    for pt in mesh.points:
        dist = lineseg_dist(pt, base, end)
        if dist < tolerance:
            pts.append(pt)

    slices = []
    for slice_pt in pts:
        slice = volume.slice(origin=slice_pt, normal=ave_direction)
        if slice.number_of_points > 0:
            slices.append(slice)

    ave_slice = slices[0].copy()

    for arrayname in ave_slice.array_names:
        ave_slice[arrayname] = ave_slice[arrayname] * 0

    for sl in slices:
        for arrayname in sl.array_names:
            ave_slice[arrayname] += sl[arrayname]

    for arrayname in ave_slice.array_names:
        ave_slice[arrayname] = ave_slice[arrayname] * 1 / len(slices)

    ave_slice = ave_slice.cell_data_to_point_data()

    if verbose:
        p = pv.Plotter()
        p.add_mesh(pv.PolyData(np.array(pts)))
        p.add_mesh(volume, show_edges=True, opacity=0.1)
        for sl in slices:
            if sl.number_of_cells > 0:
                p.add_mesh(sl, opacity=0.1)
        p.show()

    return ave_slice


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
