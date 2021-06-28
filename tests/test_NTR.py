#!/usr/bin/env python

"""Tests for `NTR` package."""

import os
import numpy as np
import pyvista as pv

from NTR.utils.functions import yaml_dict_read
from NTR.utils.geom_functions import sortProfilePoints
from NTR.utils.geom_functions import sortProfilePoints_meshing
from NTR.utils.geom_functions import calcConcaveHull


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

    xs, ys = calcConcaveHull(xs_raw, ys_raw, 1)

    assert len(xs) == len(xs_raw)
    assert any([xi in xs_raw for xi in xs])
    assert any([yi in ys_raw for yi in ys])

    polygon = pv.Polygon()
    polygon.rotate_z(np.random.randint(0,360))
    polyedges = polygon.extract_feature_edges()
    polypoints = polyedges.points
    np.random.shuffle(polypoints)
    xs_raw = polypoints[:, 0]
    ys_raw = polypoints[:, 1]

    xs, ys = calcConcaveHull(xs_raw, ys_raw, 1)

    assert len(xs) == len(xs_raw)
    assert any([xi in xs_raw for xi in xs])
    assert any([yi in ys_raw for yi in ys])

def test_profilePoints():
    print("test")
    coords_df = np.loadtxt(os.path.join("testdata_profilePoints", "profile_pointcloud.txt"))

    x_raw = coords_df[:, 0]
    y_raw = coords_df[:, 1]

    sortProfilePoints_meshing(x_raw, y_raw, 0.1)
    # sortProfilePoints(x_raw, y_raw, alpha=0.01)
