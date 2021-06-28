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


def test_concaveHull_box():
    square = pv.Plane()
    square.rotate_z(45)
    boxedges = square.extract_feature_edges()
    boxpoints = boxedges.points

    xs_raw = boxpoints[:, 0]
    ys_raw = boxpoints[:, 1]

    xs, ys = calcConcaveHull(xs_raw, ys_raw, 1)
    xs_soll = [-0.56568545, -0.49497476, -0.42426407, -0.35355338, -0.28284273, -0.21213204, -0.14142135,
               -0.070710674, 0.0, 0.070710674, 0.14142135, 0.21213204, 0.28284273, 0.35355338, 0.42426407,
               0.49497476, 0.56568545, 0.6363961, 0.70710677, 0.6363961, 0.56568545, 0.49497476, 0.42426407,
               0.35355338, 0.28284273, 0.21213204, 0.14142135, 0.070710674, 0.0, -0.070710674, -0.14142135,
               -0.21213204, -0.28284273, -0.35355338, -0.42426407, -0.49497476, -0.56568545, -0.6363961, -0.70710677,
               -0.6363961]
    ys_soll = [-0.14142135, -0.21213204, -0.28284273, -0.35355338, -0.42426407, -0.49497476, -0.56568545, -0.6363961,
               -0.70710677, -0.6363961, -0.56568545, -0.49497476, -0.42426407, -0.35355338, -0.28284273, -0.21213204,
               -0.14142135, -0.070710674, 0.0, 0.070710674, 0.14142135, 0.21213204, 0.28284273, 0.35355338,
               0.42426407, 0.49497476, 0.56568545, 0.6363961, 0.70710677, 0.6363961, 0.56568545, 0.49497476,
               0.42426407, 0.35355338, 0.28284273, 0.21213204, 0.14142135, 0.070710674, 0.0, -0.070710674]

    assert not any(np.isclose(ys, ys_soll) == False)
    assert not any(np.isclose(xs, xs_soll) == False)


def test_profilePoints():
    print("test")
    coords_df = np.loadtxt(os.path.join("testdata_profilePoints", "profile_pointcloud.txt"))

    x_raw = coords_df[:, 0]
    y_raw = coords_df[:, 1]

    sortProfilePoints_meshing(x_raw, y_raw, 0.1)
    # sortProfilePoints(x_raw, y_raw, alpha=0.01)
