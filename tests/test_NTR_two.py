import os
import pytest
import warnings

import NTR
from NTR.preprocessing.create_geom import run_create_geometry
from NTR.utils.functions import run_igg_meshfuncs, yaml_dict_read
from NTR.preprocessing.create_simcase import create_simulationcase, create_parastudsims

ntrpath = os.path.abspath(os.path.dirname(NTR.__file__))

base = os.path.join(ntrpath, "..", "examples")

examples_compressor_kurth = os.path.join(base, "CascadeCase_compressor_kurth", "case_settings.yml")
examples_gwkles = os.path.join(base, "CascadeCase_gwk_les", "case_settings.yml")
examples_gwkras = os.path.join(base, "CascadeCase_gwk_rans", "case_settings.yml")
examples_naca6510 = os.path.join(base, "CascadeCase_gwk_naca6510", "case_settings.yml")
examples_nacagen = os.path.join(base, "CascadeCase_NACA_airfoilgenerator", "case_settings.yml")
examples_turbine_seehausen = os.path.join(base, "CascadeCase_turbine_seehausen", "case_settings.yml")
examples_gwkras_trace = os.path.join(base, "CascadeCase_gwk_rans_trace", "case_settings.yml")
examples_gwkras_trace_parastud = os.path.join(base, "CascadeCase_gwk_rans_trace_parastud", "case_settings.yml")


igg_settings = yaml_dict_read(os.path.join(ntrpath, "utils", "externals", "externals_settings.yml"))
iggdir = igg_settings["igg"]["install_directory"]
igginstance = igg_settings["igg"]["executable"]
iggpath = os.path.join(iggdir, igginstance)

ON_CI = 'CI' in os.environ

def run_mesh_wrapper(case_yml):
    if os.path.isfile(iggpath):
        run_igg_meshfuncs(case_yml)
    else:
        warnings.warn("igg-path not set correctly or igg is not installed. this test will be skipped")

@pytest.mark.skipif(ON_CI, reason="do not run in continuous integration environment due to limited ressources")
def test_example_compressor_kurth():
    run_create_geometry(examples_compressor_kurth)
    run_mesh_wrapper(examples_compressor_kurth)


@pytest.mark.skipif(ON_CI, reason="do not run in continuous integration environment due to limited ressources")
def test_example_gwkles():
    run_create_geometry(examples_gwkles)
    run_mesh_wrapper(examples_gwkles)
    create_simulationcase(examples_gwkles)


@pytest.mark.skipif(ON_CI, reason="do not run in continuous integration environment due to limited ressources")
def test_example_gwkras():
    run_create_geometry(examples_gwkras)
    run_mesh_wrapper(examples_gwkras)
    create_simulationcase(examples_gwkras)


@pytest.mark.skipif(ON_CI, reason="do not run in continuous integration environment due to limited ressources")
def test_example_naca6510():
    run_create_geometry(examples_naca6510)
    run_mesh_wrapper(examples_naca6510)


@pytest.mark.skipif(ON_CI, reason="do not run in continuous integration environment due to limited ressources")
def test_example_nacagen():
    run_create_geometry(examples_nacagen)
    run_mesh_wrapper(examples_nacagen)


@pytest.mark.skipif(ON_CI, reason="do not run in continuous integration environment due to limited ressources")
def test_example_turbine_seehausen():
    run_create_geometry(examples_turbine_seehausen)
    run_mesh_wrapper(examples_turbine_seehausen)


@pytest.mark.skipif(ON_CI, reason="do not run in continuous integration environment due to limited ressources")
def test_example_gwkras_trace():
    run_create_geometry(examples_gwkras_trace)
    run_mesh_wrapper(examples_gwkras_trace)
    create_simulationcase(examples_gwkras_trace)


@pytest.mark.skipif(ON_CI, reason="do not run in continuous integration environment due to limited ressources")
def test_example_gwkras_trace_parastud():
    run_create_geometry(examples_gwkras_trace_parastud)
    run_mesh_wrapper(examples_gwkras_trace_parastud)
    create_parastudsims(examples_gwkras_trace_parastud)
