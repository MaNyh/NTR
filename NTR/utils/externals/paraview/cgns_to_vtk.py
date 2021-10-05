import os
import pyvista as pv
import shutil

import NTR
from NTR.utils.filehandling import yaml_dict_read

externals = yaml_dict_read(os.path.join(os.path.dirname(NTR.__file__),"utils","externals","externals_settings.yml"))


cgns_to_read = os.path.join("D://","01_CascadeCase_gwk_ras_trace_initialfit","02_Simcase","output","cgns","TRACE.cgns")
target = os.path.join("D://","01_CascadeCase_gwk_ras_trace_initialfit","03_Solution","solution.vtk")


def convert_cgns_to_vtk(cgns_to_read,target):
    tmpcgns = "tmp.cgns"
    scriptname = "paraview_template_cgns_to_vtk.py"
    shutil.copyfile(cgns_to_read, tmpcgns)
    os.system(externals["paraview"]["pvpython"] + " " + scriptname)
    os.remove(tmpcgns)
    shutil.copyfile("solution.vtk", target)
    os.remove("solution.vtk")


convert_cgns_to_vtk(cgns_to_read,target)
