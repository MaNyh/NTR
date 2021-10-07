import os
import shutil

import NTR
from NTR.utils.filehandling import yaml_dict_read

externals = yaml_dict_read(os.path.join(os.path.dirname(NTR.__file__),"utils","externals","externals_settings.yml"))


def convert_cgns_to_vtk(cgns_to_read, target):

    tmpcgns = "tmp.cgns"
    tmpvtk = "solution.vtk"
    #tmpcgnspath = os.path.join(cwd,tmpcgns)
    scriptname = "paraview_template_cgns_to_vtk.py"
    #if os.path.isfile(os.path.join(cwd,))
    shutil.copyfile(cgns_to_read, os.path.join(target,tmpcgns))
    os.system(externals["paraview"]["pvpython"] + " " + os.path.join(os.path.dirname(__file__),scriptname))
    shutil.copyfile(os.path.join(os.path.dirname(__file__),"solution.vtk"), os.path.join(target,tmpvtk))
    os.remove(os.path.join(target,tmpcgns))
