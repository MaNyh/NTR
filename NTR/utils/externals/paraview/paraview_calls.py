import os
import shutil
import tempfile

import NTR
from NTR.database.case_dirstructure import casedirs
from NTR.utils.filehandling import yaml_dict_read

externals = yaml_dict_read(os.path.join(os.path.dirname(NTR.__file__),"utils","externals","externals_settings.yml"))

def extract_patches(mesh,inlet,outlet,blade):
    temp_dir = tempfile.TemporaryDirectory()
    scriptname = "traceblade.py"
    scriptsource = os.path.join(os.path.dirname(__file__),scriptname)
    scripttarget = os.path.join(temp_dir.name,scriptname)

    with open(scriptsource, "rt") as fin:
        with open(scripttarget, "w") as fout:
            for line in fin:
                nline = line.replace("<var MESH var>", os.path.abspath(mesh))
                nline = nline.replace("<var INLET var>", os.path.abspath(inlet))
                nline = nline.replace("<var OUTLET var>", os.path.abspath(outlet))
                nline = nline.replace("<var BLADE var>", os.path.abspath(blade))
                fout.write(nline)
    os.system(externals["paraview"]["pvpython"] + " " + scripttarget)
    temp_dir.cleanup()


def convert_cgns_to_vtk(cgns_to_read, target):

    tmpcgns = "tmp.cgns"
    tmpvtk = "solution.vtk"
    #tmpcgnspath = os.path.join(cwd,tmpcgns)
    scriptname = "paraview_template_cgns_to_vtk.py"
    shutil.copyfile(cgns_to_read, os.path.join(os.path.dirname(__file__),tmpcgns))
    os.system(externals["paraview"]["pvpython"] + " " + os.path.join(os.path.dirname(__file__),scriptname))
    shutil.copyfile(os.path.join(os.path.dirname(__file__),"solution.vtk"), target)
    os.remove(os.path.join(target,os.path.join(os.path.dirname(__file__),tmpcgns)))
    os.remove(os.path.join(target,os.path.join(os.path.dirname(__file__),tmpvtk)))


def paraview_convert_cgns_to_vtk(settings_yml):
    settings = yaml_dict_read(settings_yml)
    case_path = os.path.dirname(settings_yml)

    assert "post_settings" in settings.keys(), "post_settings missing"
    assert "paraview_convert_cgns_to_vtk" in settings["post_settings"].keys(), "paraview_convert_cgns_to_vtk settings missing"
    assert "cgns" in settings["post_settings"]["paraview_convert_cgns_to_vtk"].keys(), "no cgns to convert defined in settings"
    path_to_volmesh = os.path.join(case_path, settings["post_settings"]["paraview_convert_cgns_to_vtk"]["cgns"])

    assert os.path.isfile(path_to_volmesh), "file " + path_to_volmesh + " does not exist"

    target_name = settings["post_settings"]["paraview_convert_cgns_to_vtk"]["vtk_name"]
    convert_cgns_to_vtk(path_to_volmesh, os.path.join(case_path, casedirs["solution"], target_name))
    return 0
