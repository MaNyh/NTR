import os

import NTR
from NTR.preprocessing.create_geom import run_create_geometry
from NTR.utils.functions import run_igg_meshfuncs, yaml_dict_read
from NTR.preprocessing.create_simcase import create_simulationcase

ntrpath = os.path.abspath(os.path.dirname(NTR.__file__))
examples_compressor_kurth = os.path.join(ntrpath,"..", "examples", "CascadeCase_compressor_kurth", "case_settings.yml")
examples_gwkles = os.path.join(ntrpath,"..", "examples", "CascadeCase_gwk_les", "case_settings.yml")
examples_gwkras = os.path.join(ntrpath,"..", "examples", "CascadeCase_gwk_rans", "case_settings.yml")
examples_naca6510 = os.path.join(ntrpath,"..", "examples", "CascadeCase_gwk_naca6510", "case_settings.yml")
examples_nacagen = os.path.join(ntrpath,"..", "examples", "CascadeCase_NACA_airfoilgenerator", "case_settings.yml")
examples_turbine_seehausen = os.path.join(ntrpath,"..", "examples", "CascadeCase_turbine_seehausen", "case_settings.yml")

igg_settings = yaml_dict_read("../NTR/utils/externals/externals_settings.yml")
iggdir = igg_settings["igg"]["install_directory"]
igginstance = igg_settings["igg"]["executable"]
iggpath = os.path.join(iggdir, igginstance)

@pytest.mark.skipif(ON_CI = 'CI' in os.environ, reason="do not run this on a server-venv")
def test_example_compressor_kurth():
    run_create_geometry(examples_compressor_kurth)
    if os.path.isfile(iggpath):
        run_igg_meshfuncs(examples_compressor_kurth)
    else:
        print("igg-path is not right")

@pytest.mark.skipif(ON_CI = 'CI' in os.environ, reason="do not run this on a server-venv")
def test_example_gwkles():
    if os.path.isfile(iggpath):
        run_igg_meshfuncs(examples_gwkles)
    else:
        print("igg-path is not right")
    create_simulationcase(examples_gwkles)

@pytest.mark.skipif(ON_CI = 'CI' in os.environ, reason="do not run this on a server-venv")
def test_example_gwkras():
    run_create_geometry(examples_gwkras)

    if os.path.isfile(iggpath):
        run_igg_meshfuncs(examples_gwkras)
    else:
        print("igg-path is not right")
    create_simulationcase(examples_gwkras)

@pytest.mark.skipif(ON_CI = 'CI' in os.environ, reason="do not run this on a server-venv")
def test_example_naca6510():
    run_create_geometry(examples_naca6510)

    if os.path.isfile(iggpath):
        run_igg_meshfuncs(examples_naca6510)
    else:
        print("igg-path is not right")

@pytest.mark.skipif(ON_CI = 'CI' in os.environ, reason="do not run this on a server-venv")
def test_example_nacagen():
    run_create_geometry(examples_nacagen)

    if os.path.isfile(iggpath):
        run_igg_meshfuncs(examples_nacagen)
    else:
        print("igg-path is not right")

@pytest.mark.skipif(ON_CI = 'CI' in os.environ, reason="do not run this on a server-venv")
def test_example_turbine_seehausen():
    run_create_geometry(examples_turbine_seehausen)

    if os.path.isfile(iggpath):
        run_igg_meshfuncs(examples_turbine_seehausen)
    else:
        print("igg-path is not right")
