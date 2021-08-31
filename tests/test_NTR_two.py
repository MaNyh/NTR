import os


from NTR.preprocessing.create_geom import run_create_geometry
from NTR.utils.functions import run_igg_meshfuncs
from NTR.preprocessing.create_simcase import create_simulationcase

examples_compressor_kurth = os.path.join("..", "examples", "CascadeCase_compressor_kurth", "case_settings.yml")
examples_gwkles = os.path.join("..", "examples", "CascadeCase_gwk_les", "case_settings.yml")
examples_gwkras = os.path.join("..", "examples", "CascadeCase_gwk_rans", "case_settings.yml")
examples_naca6510 = os.path.join("..", "examples", "CascadeCase_gwk_naca6510", "case_settings.yml")
examples_nacagen = os.path.join("..", "examples", "CascadeCase_NACA_airfoilgenerator", "case_settings.yml")
examples_turbine_seehausen = os.path.join("..", "examples", "CascadeCase_turbine_seehausen", "case_settings.yml")


def test_example_compressor_kurth():
    #run_create_geometry(examples_compressor_kurth)
    run_igg_meshfuncs(examples_compressor_kurth)


def test_example_gwkles():
    run_create_geometry(examples_gwkles)
    run_igg_meshfuncs(examples_gwkles)
    create_simulationcase(examples_gwkles)


def test_example_gwkras():
    run_create_geometry(examples_gwkras)
    run_igg_meshfuncs(examples_gwkras)
    create_simulationcase(examples_gwkras)


def test_example_naca6510():
    run_create_geometry(examples_naca6510)
    run_igg_meshfuncs(examples_naca6510)


def test_example_nacagen():
    run_create_geometry(examples_nacagen)
    run_igg_meshfuncs(examples_nacagen)


def test_example_turbine_seehausen():
    run_create_geometry(examples_turbine_seehausen)
    run_igg_meshfuncs(examples_turbine_seehausen)
