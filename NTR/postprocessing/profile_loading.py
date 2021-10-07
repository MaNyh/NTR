import os
import tempfile

from NTR.utils.externals.paraview.cgns_to_vtk import convert_cgns_to_vtk
from NTR.utils.geom_functions.pyvista_utils import load_mesh

numerical_mesh = "D:\\01_CascadeCase_gwk_ras_trace_initialfit\\02_Simcase\\output\\cgns\\TRACE.cgns"

def calc_loading_volmesh(path_to_volmesh):
    assert os.path.isfile(path_to_volmesh),"file "+ path_to_volmesh + "does not exist"
    temp_dir = tempfile.TemporaryDirectory()

    mesh_ext = os.path.splitext(path_to_volmesh)[-1]
    volmesh = None

    if mesh_ext == ".cgns":
        tmpmesh = "solution.vtk"
        convert_cgns_to_vtk(path_to_volmesh,temp_dir.name)
        volmesh = load_mesh(os.path.join(temp_dir.name,tmpmesh))
        temp_dir.cleanup()
    elif mesh_ext == ".vtk":
        volmesh(path_to_volmesh)

    return volmesh

volmesh = calc_loading_volmesh(numerical_mesh)
