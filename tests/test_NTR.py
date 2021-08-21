#!/usr/bin/env python

"""Tests for `NTR` package."""

import numpy as np
import pyvista as pv

from NTR.utils.filehandling import yaml_dict_read
from NTR.utils.geom_functions import sortProfilePoints
from NTR.utils.geom_functions import calcConcaveHull
from NTR.utils.geom_functions import extract_vk_hk


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


def test_extract_vk_hk():
    """
    tests a NACA0009 symmetric profile in a random angle as a minimal example.
    hk is point id 9, vk is id 34
    :return:
    """
    points2d = np.array([[1.00000, 0.0],
                         [0.99572, 0.00057],
                         [0.98296, 0.00218],
                         [0.96194, 0.00463],
                         [0.93301, 0.00770],
                         [0.89668, 0.01127],
                         [0.85355, 0.01522],
                         [0.80438, 0.01945],
                         [0.75000, 0.02384],
                         [0.69134, 0.02823],
                         [0.62941, 0.03247],
                         [0.56526, 0.03638],
                         [0.50000, 0.03978],
                         [0.43474, 0.04248],
                         [0.37059, 0.04431],
                         [0.33928, 0.04484],
                         [0.30866, 0.04509],
                         [0.27886, 0.04504],
                         [0.25000, 0.04466],
                         [0.22221, 0.04397],
                         [0.19562, 0.04295],
                         [0.17033, 0.04161],
                         [0.14645, 0.03994],
                         [0.12408, 0.03795],
                         [0.10332, 0.03564],
                         [0.08427, 0.03305],
                         [0.06699, 0.03023],
                         [0.05156, 0.02720],
                         [0.03806, 0.02395],
                         [0.02653, 0.02039],
                         [0.01704, 0.01646],
                         [0.00961, 0.01214],
                         [0.00428, 0.00767],
                         [0.00107, 0.00349],
                         [0.0, 0.0],
                         [0.00107, -0.00349],
                         [0.00428, -0.00767],
                         [0.00961, -0.01214],
                         [0.01704, -0.01646],
                         [0.02653, -0.02039],
                         [0.03806, -0.02395],
                         [0.05156, -0.02720],
                         [0.06699, -0.03023],
                         [0.08427, -0.03305],
                         [0.10332, -0.03564],
                         [0.12408, -0.03795],
                         [0.14645, -0.03994],
                         [0.17033, -0.04161],
                         [0.19562, -0.04295],
                         [0.22221, -0.04397],
                         [0.25000, -0.04466],
                         [0.27886, -0.04504],
                         [0.30866, -0.04509],
                         [0.33928, -0.04484],
                         [0.37059, -0.04431],
                         [0.43474, -0.04248],
                         [0.50000, -0.03978],
                         [0.56526, -0.03638],
                         [0.62941, -0.03247],
                         [0.69134, -0.02823],
                         [0.75000, -0.02384],
                         [0.80438, -0.01945],
                         [0.85355, -0.01522],
                         [0.89668, -0.01127],
                         [0.93301, -0.00770],
                         [0.96194, -0.00463],
                         [0.98296, -0.00218],
                         [0.99572, -0.00057]])

    points = np.stack((points2d[:, 0], points2d[:, 1], np.zeros(len(points2d)))).T

    profilepoints = pv.PolyData(points)

    random_angle = np.random.randint(-40, 40)
    profilepoints.rotate_z(random_angle)

    origPoly = pv.PolyData(profilepoints)
    sortedPoly = pv.PolyData(profilepoints)
    ind_hk, ind_vk, veronoi_mid = extract_vk_hk(origPoly, sortedPoly, verbose=False)
    assert ind_hk == 0, "wrong hk-index found"
    assert ind_vk == 34, "wrong vk-index found"
