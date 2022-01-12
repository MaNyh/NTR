#!/usr/bin/env python

"""Tests for `NTR` package."""
import os
import shutil

import numpy as np
import pyvista as pv
import yaml

import NTR
from NTR.preprocessing.create_simcase import create_simulationcase, find_vars_opts, read_parastudyaml, \
    get_common_association, copy_template, create_simdirstructure
from NTR.utils.dicthandling import nested_dict_pairs_iterator, merge
from NTR.utils.filehandling import yaml_dict_read, write_yaml_dict, write_pickle, get_directory_structure
from NTR.utils.geom_functions.pointcloud import calcConcaveHull
from NTR.utils.geom_functions.profileparas import extract_vk_hk, sortProfilePoints, extractSidePolys, midline_from_sides
from NTR.utils.geom_functions.spline import splineCurvature
from NTR.database.case_dirstructure import casedirs
from NTR.postprocessing.generic.spatial_average import vol_to_plane, vol_to_line
from NTR.postprocessing.generic.integralscales_from_signal import integralscales_from_timeseries


def test_yamlDictRead(tmpdir):
    """
    tmpdir ist eine spezialvariable die durch pytest erkannt wird (ist ein PYTHONPATH-objekt)
    """
    test_file = tmpdir / "test.yaml"
    with open(test_file, "w") as handle:
        handle.write("test_key: True\n")
    assert yaml_dict_read(test_file) == {"test_key": True}


def test_splineCurvature():
    radius = np.random.rand()
    curve = 1 / radius
    circle = pv.CircularArc((-radius, 0, 0), (radius, 0, 0), (0, 0, 0))
    circle.rotate_z(-180,inplace=False)
    curvature_return = splineCurvature(circle.points[::, 0], circle.points[::, 1])
    assert not any(abs(curvature_return[3:-3] - curve) > 0.05), "spline curvature calculation failed"


def test_calcConcaveHull():
    """
    in these simple geometries, each point must be found by calcConcaveHull
    """
    square = pv.Plane()
    boxedges = square.extract_feature_edges()

    boxedges.rotate_z(np.random.randint(0, 360),inplace=False)
    boxpoints = boxedges.points

    np.random.shuffle(boxpoints)

    xs_raw = boxpoints[:, 0]
    ys_raw = boxpoints[:, 1]

    xs, ys = calcConcaveHull(xs_raw, ys_raw, 10)

    assert len(xs) == len(xs_raw)
    assert any([yi in ys_raw for yi in ys])

    polygon = pv.Polygon()
    polygon.rotate_z(np.random.randint(0, 360),inplace=False)
    polyedges = polygon.extract_feature_edges()
    polypoints = polyedges.points
    np.random.shuffle(polypoints)
    xs_raw = polypoints[:, 0]
    ys_raw = polypoints[:, 1]

    xs, ys = calcConcaveHull(xs_raw, ys_raw, 10)

    assert len(xs) == len(xs_raw)
    assert any([yi in ys_raw for yi in ys])


def test_profilePoints():
    ellipse = pv.ParametricEllipsoid(1, np.random.rand(), np.random.rand())
    ellipse = ellipse.slice(normal="z", origin=(0, 0, 0))
    test_area = pv.PolyData(ellipse.points)
    test_area = test_area.delaunay_2d()
    test_area = test_area.compute_cell_sizes()

    tarea = sum(test_area["Area"])

    ellipse.rotate_z(np.random.randint(0, 360),inplace=False)
    shape = ellipse.extract_geometry()
    pts = shape.points
    xs = pts[:, 0]
    ys = pts[:, 1]

    x_ss, y_ss, x_ps, y_ps = sortProfilePoints(xs, ys, 50)

    pts_x = x_ss + x_ps
    pts_y = y_ss + y_ps

    reconstructed_poly = pv.PolyData(np.stack((pts_x, pts_y, np.zeros(len(pts_x)))).T)
    reconstructed_face = reconstructed_poly.delaunay_2d()
    reconstructed_face = reconstructed_face.compute_cell_sizes()
    reconstructed_area = sum(reconstructed_face["Area"])

    assert np.isclose(reconstructed_area, tarea)


def test_extract_vk_hk(verbose=False):
    """
    tests a NACA  profile in a random angle as a minimal example.
    :return:
    """
    from NTR.database.data_generators.naca_airfoil_creator import naca

    res = 400

    # d1,d2,d3,d4 = np.random.randint(0,9),np.random.randint(0,9),np.random.randint(0,9),np.random.randint(0,9)
    # digitstring = str(d1)+str(d2)+str(d3)+str(d4)
    # manifold problems with other profiles with veronoi-mid and other unoptimized code. therefor tests only 0009
    X, Y = naca("6409", res, finite_TE=False, half_cosine_spacing=True)
    ind_hk_test = 0
    ind_vk_test = res

    points = np.stack((X[:-1], Y[:-1], np.zeros(res * 2))).T

    profilepoints = pv.PolyData(points)

    random_angle = np.random.randint(-40, 40)
    profilepoints.rotate_z(random_angle,inplace=False)

    sortedPoly = pv.PolyData(profilepoints)
    ind_hk, ind_vk = extract_vk_hk(sortedPoly, verbose=verbose)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(sortedPoly.points[ind_hk], color="yellow", point_size=20)
        p.add_mesh(sortedPoly.points[ind_vk], color="red", point_size=20)
        p.add_mesh(sortedPoly)
        p.show()

    assert ind_hk == ind_hk_test, "wrong hk-index chosen"
    assert ind_vk == ind_vk_test, "wrong vk-index chosen"


def test_extractSidePolys(verbose=False):
    from NTR.database.data_generators.naca_airfoil_creator import naca
    d1, d2, d3, d4 = np.random.randint(0, 9), np.random.randint(0, 9), np.random.randint(0, 9), np.random.randint(0, 9)
    digitstring = str(d1) + str(d2) + str(d3) + str(d4)

    res = 240
    X, Y = naca(digitstring, res, finite_TE=False, half_cosine_spacing=True)
    ind_hk = 0
    ind_vk = res
    points = np.stack((X[:-1], Y[:-1], np.zeros(res * 2) - 1)).T

    poly = pv.PolyData(points)
    ssPoly, psPoly = extractSidePolys(ind_hk, ind_vk, poly, verbose)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(poly.points[ind_hk], color="yellow", point_size=20)
        p.add_mesh(poly.points[ind_vk], color="yellow", point_size=20)
        p.add_mesh(ssPoly, color="blue")
        p.add_mesh(psPoly, color="red")
        p.show()

    assert ssPoly.number_of_points == psPoly.number_of_points, "number of sidepoints are not equal"


def test_midline_from_sides(verbose=False):
    from NTR.utils.mathfunctions import vecAbs
    from NTR.database.data_generators.naca_airfoil_creator import naca

    res = 240
    X, Y = naca('0009', res, finite_TE=False, half_cosine_spacing=True)
    ind_hk = 0
    ind_vk = res

    points = np.stack((X[:-1], Y[:-1], np.zeros(res * 2) - 1)).T
    poly = pv.PolyData(points)
    ssPoly, psPoly = extractSidePolys(ind_hk, ind_vk, poly)

    mids = midline_from_sides(ind_hk, ind_vk, poly.points, psPoly, ssPoly)

    length = mids.length
    testlength = vecAbs(ssPoly.points[0] - ssPoly.points[-1])

    if verbose:
        p = pv.Plotter()
        p.add_mesh(mids, color="yellow", point_size=20)
        p.add_mesh(ssPoly, color="blue")
        p.add_mesh(psPoly, color="red")
        p.show()

    assert length == testlength, "midline not accurate"


def test_create_simulationcase(tmpdir):
    ntrpath = os.path.abspath(os.path.dirname(NTR.__file__))
    case_templates = os.listdir(os.path.join(ntrpath, "database", "case_templates"))

    case_structure_templates = {}
    templates_basedir = os.path.join(ntrpath, "database", "case_templates")

    for cname in case_templates:
        cstruct = get_directory_structure(os.path.join(templates_basedir, cname))
        case_structure_templates[cname] = cstruct[cname]

    for indx in range(len(list(case_structure_templates.keys()))):
        case_type = list(case_structure_templates.keys())[indx]
        case_structure_template = case_structure_templates[case_type]

        common_dir = get_common_association(case_type)
        create_simdirstructure(case_structure_template, tmpdir)
        copy_template(case_type, case_structure_template, tmpdir)
        all_pairs = list(nested_dict_pairs_iterator(case_structure_template))

        if common_dir:
            common_base = os.path.join(ntrpath, "database", "common_files", common_dir)

            swap_pairs = [[idx, i] for idx, i in enumerate(all_pairs) if i[-2][-6:] == "common"]

            for idx, sw in swap_pairs:
                dirs = sw[:-2]
                file = sw[-2]
                # todo: overthink, weather not (common_base,*dirs,file) would not be a better suiting structure for these templatefiles
                source = os.path.join(common_base, file)
                target = os.path.join(tmpdir, *dirs, file.replace(".common", ""))
                all_pairs_helper = [list(i) for i in all_pairs]
                all_pairs_helper[idx][-2] = list(all_pairs_helper[idx])[-2].replace(".common", "")
                all_pairs = all_pairs_helper
                shutil.copyfile(source, target)
                os.remove(os.path.join(tmpdir, *dirs, file))

        case_structure_var = find_vars_opts(case_structure_template, "var", all_pairs, tmpdir)
        case_structure_opt = find_vars_opts(case_structure_template, "opt", all_pairs, tmpdir)
        case_structure = merge(case_structure_var, case_structure_opt)
        case_structlist = list(nested_dict_pairs_iterator(case_structure))
        variables = [i[-2] for i in list(case_structlist) if i[-1] == "var"]
        options = [i[-2] for i in list(case_structlist) if i[-1] == "opt"]

        test_dict = {"case_settings": {"case_type": case_type,
                                       "name": "testcase",
                                       "type": "simulation",
                                       "job": {
                                           "job_script": "submit_sp_pbs_rrzn",
                                           "job_nodes": "1",
                                           "job_ppn": "1",
                                           "job_mail": "asd@asd.asd",
                                           "job_mem": "12",
                                       },
                                       "description": "this is a test-function",
                                       "sub_cmd": "qsub", },
                     "simcase_settings": {"variables": {},
                                          "options": {}},
                     "simcase_optiondef": {},
                     }

        for var in variables:
            test_dict["simcase_settings"]["variables"][var] = "1"
        for opt in options:
            test_dict["simcase_settings"]["options"][opt] = "1"
            test_dict["simcase_optiondef"][opt] = "1"

        test_file = tmpdir / "test_create_simulationcase.yml"
        test_geo_dict = {}
        if not os.path.isdir(tmpdir / casedirs["data"]):
            os.mkdir(tmpdir / casedirs["data"])
        test_geo_file = tmpdir / os.path.join(casedirs["data"], "geometry.pkl")
        write_yaml_dict(test_file, test_dict)
        write_pickle(test_geo_file, test_geo_dict)

        create_simulationcase(test_file)
        for root, dirs, files in os.walk(tmpdir):
            for f in files:
                os.unlink(os.path.join(root, f))
            for d in dirs:
                shutil.rmtree(os.path.join(root, d))


def test_read_fullfactorparastud_yaml(tmpdir):
    keyword = "testsettings"

    listname = "list"
    listvalues = [0, 1]

    integername = "integer"
    integervalue = 1

    test_dict = {keyword: {listname: listvalues,
                           integername: integervalue}
                 }

    test_file = tmpdir / "test.yaml"
    with open(test_file, "w") as handle:
        yaml.dump(test_dict, handle, default_flow_style=False)

    testsets = read_parastudyaml(test_file)
    sets = [{keyword: {integername: integervalue, listname: listvalues[0]}},
            {keyword: {integername: integervalue, listname: listvalues[1]}}]

    assert testsets == sets, "parametrized dict is not interpreted right"


def test_vol_to_plane():
    dirs = ["x", "y", "z"]
    dir_idx_dict = np.random.choice(range(len(dirs)))
    dir = dirs[dir_idx_dict]

    values = np.linspace(0, 10, 1000).reshape((20, 5, 10))
    values.shape

    # Create the spatial reference
    grid = pv.UniformGrid()

    # Set the grid dimensions: shape + 1 because we want to inject our values on
    #   the CELL data
    grid.dimensions = np.array(values.shape) + 1

    # Edit the spatial reference
    grid.origin = (100, 33, 55.6)  # The bottom left corner of the data set
    grid.spacing = (1, 5, 2)  # These are the cell sizes along each axis
    grid.clear_arrays()
    grid.point_arrays["var"] = grid.points[::, dir_idx_dict]
    grid = grid.point_data_to_cell_data()
    plane = vol_to_plane(grid, dir)

    grid_cl_high = grid.bounds[dir_idx_dict * 2 + 1]
    grid_cl_low = grid.bounds[dir_idx_dict * 2]

    meanval = (grid_cl_high + grid_cl_low) / 2
    assert len(np.where(np.isclose(plane["var"], meanval))[0]) == plane.number_of_points


def test_vol_to_line():
    dirs = ["x", "y", "z"]
    dir_idx_dict = np.random.choice(range(len(dirs)))
    dir = dirs[dir_idx_dict]

    values = np.linspace(0, 10, 1000).reshape((20, 5, 10))
    values.shape

    # Create the spatial reference
    grid = pv.UniformGrid()

    # Set the grid dimensions: shape + 1 because we want to inject our values on
    #   the CELL data
    grid.dimensions = np.array(values.shape) + 1

    # Edit the spatial reference
    grid.origin = (100, 33, 55.6)  # The bottom left corner of the data set
    grid.spacing = (1, 5, 2)  # These are the cell sizes along each axis
    grid.clear_arrays()
    grid.point_arrays["var"] = grid.points[::, dir_idx_dict]
    grid = grid.point_data_to_cell_data()
    pts, var = vol_to_line(grid, dir)

    grid_cl_high = grid.bounds[dir_idx_dict * 2 + 1]
    grid_cl_low = grid.bounds[dir_idx_dict * 2]

    meanval = (grid_cl_high + grid_cl_low) / 2
    meanvar = np.mean(var["var"])
    assert np.isclose(meanvar, meanval)


def test_integralscales():
    time = np.arange(0, 1000, 0.1)
    fluctationstrength = 1
    meanvalue = 3
    frequency = 2
    signal = np.sin(time * frequency) * fluctationstrength + meanvalue
    mean = np.mean(signal)
    fluctation = signal - mean
    timescale, lenghtscale = integralscales_from_timeseries(mean, fluctation, time)
    assert np.isclose(timescale, frequency ** -1,
                      rtol=0.075), "calculated timescale result out of tolerance. something's wrong"
    assert np.isclose(lenghtscale, timescale * meanvalue,
                      rtol=0.075), "calculated length scale out of tolerance. something's wrong"
    return 0
