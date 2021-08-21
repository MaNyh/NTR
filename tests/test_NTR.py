#!/usr/bin/env python

"""Tests for `NTR` package."""

import numpy as np
import pyvista as pv

from NTR.utils.filehandling import yaml_dict_read
from NTR.utils.geom_functions.pointcloud import calcConcaveHull
from NTR.utils.geom_functions.profileparas import extract_vk_hk, sortProfilePoints, extractSidePolys, midline_from_sides


def test_yamlDictRead(tmpdir):
    """
    tmpdir ist eine spezialvariable die durch pytest erkannt wird (ist ein PYTHONPATH-objekt)
    """
    test_file = tmpdir / "test.yaml"
    with open(test_file, "w") as handle:
        handle.write("test_key: True\n")
    assert yaml_dict_read(test_file) == {"test_key": True}


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
    tests a NACA0009 symmetric profile in a random angle as a minimal example.
    :return:
    """
    from tests.datasets_test import naca0009profile
    points2d = naca0009profile.points

    ind_hk_test = naca0009profile.ind_hk
    ind_vk_test = naca0009profile.ind_vk

    points = np.stack((points2d[:, 0], points2d[:, 1], np.zeros(len(points2d)))).T

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

    assert ind_hk == ind_hk_test, "wrong hk-index found"
    assert ind_vk == ind_vk_test, "wrong vk-index found"


def test_extractSidePolys(verbose=False):
    from tests.datasets_test import naca0009profile
    points2d = naca0009profile.points
    ind_hk = naca0009profile.ind_hk
    ind_vk = naca0009profile.ind_vk

    poly = pv.PolyData(np.stack((points2d[:, 0], points2d[:, 1], np.zeros(len(points2d)))).T)
    ssPoly, psPoly = extractSidePolys(ind_hk, ind_vk, poly)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(poly.points[ind_hk], color="yellow", point_size=20)
        p.add_mesh(poly.points[ind_vk], color="yellow", point_size=20)
        p.add_mesh(ssPoly, color="blue")
        p.add_mesh(psPoly, color="red")
        p.show()

    assert ssPoly.number_of_points == psPoly.number_of_points, "number of sidepoints are not equal, test failed"


def test_midline_from_sides(verbose=False):
    from NTR.utils.mathfunctions import vecAbs
    from tests.datasets_test import naca0009profile

    points2d = naca0009profile.points
    ind_hk = naca0009profile.ind_hk
    ind_vk = naca0009profile.ind_vk

    poly = pv.PolyData(np.stack((points2d[:, 0], points2d[:, 1], np.zeros(len(points2d)))).T)
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
