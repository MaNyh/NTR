#!/usr/bin/env python

"""Tests for `NTR` package."""
import os

import numpy as np
import pyvista as pv

import NTR
from NTR.preprocessing.create_simcase import create_simulationcase, find_vars_opts, nested_dict_pairs_iterator
from NTR.utils.filehandling import yaml_dict_read, write_yaml_dict, write_pickle, get_directory_structure
from NTR.utils.geom_functions.pointcloud import calcConcaveHull
from NTR.utils.geom_functions.profileparas import extract_vk_hk, sortProfilePoints, extractSidePolys, midline_from_sides
from NTR.utils.geom_functions.spline import splineCurvature


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
    circle.rotate_z(-180)
    curvature_return = splineCurvature(circle.points[::, 0], circle.points[::, 1])
    assert not any(abs(curvature_return[3:-3] - curve) > 0.05), "spline curvature calculation failed"


def test_calcConcaveHull():
    """
    in these simple geometries, each point must be found by calcConcaveHull
    """
    square = pv.Plane()
    boxedges = square.extract_feature_edges()

    boxedges.rotate_z(np.random.randint(0, 360))
    boxpoints = boxedges.points

    np.random.shuffle(boxpoints)

    xs_raw = boxpoints[:, 0]
    ys_raw = boxpoints[:, 1]

    xs, ys = calcConcaveHull(xs_raw, ys_raw, 10)

    assert len(xs) == len(xs_raw)
    assert any([yi in ys_raw for yi in ys])

    polygon = pv.Polygon()
    polygon.rotate_z(np.random.randint(0, 360))
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

    ellipse.rotate_z(np.random.randint(0, 360))
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

    res = 80

    # d1,d2,d3,d4 = np.random.randint(0,9),np.random.randint(0,9),np.random.randint(0,9),np.random.randint(0,9)
    # digitstring = str(d1)+str(d2)+str(d3)+str(d4)
    # manifold problems with other profiles with veronoi-mid and other unoptimized code. therefor tests only 0009
    X, Y = naca("6504", res, finite_TE=False, half_cosine_spacing=True)
    ind_hk_test = 0
    ind_vk_test = res

    points = np.stack((X[:-1], Y[:-1], np.zeros(res * 2))).T

    profilepoints = pv.PolyData(points)

    random_angle = np.random.randint(-40, 40)
    profilepoints.rotate_z(random_angle)

    origPoly = pv.PolyData(profilepoints)
    sortedPoly = pv.PolyData(profilepoints)
    ind_hk, ind_vk, veronoi_mid = extract_vk_hk(origPoly, sortedPoly, verbose=verbose)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(sortedPoly.points[ind_hk], color="red", point_size=20)
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

    case_structures = {}
    for cname in case_templates:
        cstruct = get_directory_structure(os.path.join(ntrpath, "database", "case_templates"))
        case_structures[cname] = cstruct

    for indx in range(len(list(case_structures.keys()))):
        case_type = list(case_structures.keys())[indx]
        case_structure = case_structures[case_type]["case_templates"]
        case_structure = find_vars_opts(case_structure)
        case_structlist = list(nested_dict_pairs_iterator(case_structure))
        variables = [i[-2] for i in list(case_structlist) if i[-1] == "var"]
        options = [i[-2] for i in list(case_structlist) if i[-1] == "opt"]

        test_dict = {"case_settings": {"case_type": case_type,
                                       "name": "testcase"},
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
        if not os.path.isdir(tmpdir / "04_Data"):
            os.mkdir(tmpdir / "04_Data")
        test_geo_file = tmpdir / os.path.join("04_Data", "geometry.pkl")
        write_yaml_dict(test_file, test_dict)
        write_pickle(test_geo_file, test_geo_dict)
        create_simulationcase(test_file)
