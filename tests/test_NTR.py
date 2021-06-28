#!/usr/bin/env python

"""Tests for `NTR` package."""


import pytest
from click.testing import CliRunner

from NTR import NTR
from NTR import cli

from NTR.utils.functions import yaml_dict_read



def test_yamlDictRead(tmpdir):
    """
    tmpdir ist eine spezialvariable die durch pytest erkannt wird (ist ein PYTHONPATH-objekt)
    """
    test_file = tmpdir / "test.yaml"
    with open(test_file, "w") as handle:
        handle.write("test_key: True\n")
    assert yaml_dict_read(test_file) == {"test_key": True}



