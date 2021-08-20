#!/usr/bin/env python

"""Tests for `NTR` package."""

#import os
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

    boxedges.rotate_z(np.random.randint(0,360))
    boxpoints = boxedges.points

    np.random.shuffle(boxpoints)

    xs_raw = boxpoints[:, 0]
    ys_raw = boxpoints[:, 1]

    xs, ys = calcConcaveHull(xs_raw, ys_raw, 10)

    assert len(xs) == len(xs_raw)
    assert any([yi in ys_raw for yi in ys])

    polygon = pv.Polygon()
    polygon.rotate_z(np.random.randint(0,360))
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

    points2d = np.loadtxt("symmairfoil.pts")

    points = np.stack((points2d[:,0],points2d[:,1],np.zeros(len(points2d)))).T

    profilepoints = pv.PolyData(points)

    random_angle = np.random.randint(-40,40)
    profilepoints.rotate_z(random_angle)

    origPoly = pv.PolyData(profilepoints)
    sortedPoly = pv.PolyData(profilepoints)
    ind_hk, ind_vk, veronoi_mid = extract_vk_hk(origPoly,sortedPoly,verbose=False)
    assert ind_hk == 0, "wrong hk-index found"
    assert ind_vk == 34, "wrong vk-index found"
