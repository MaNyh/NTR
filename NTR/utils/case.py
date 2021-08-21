import os

from NTR.utils.coefficients import FluidCoeffs, CascadeCoeffs
from NTR.utils.geom_functions.pyvista_utils import slice_midspan_z, load_mesh, mesh_scalar_gradients
from NTR.preprocessing.openfoam.solver_variable_dicts import solver_var_dicts


class AbstractCase:
    def __init__(self, name, vartype, simtype):
        self.name = name
        self.casedir = None
        self.mesh_dict = {}
        self.mesh_loaded_dict = {}
        self.FluidCoeffs = FluidCoeffs()

        self.sim_type = simtype
        self.var_dict = solver_var_dicts[vartype]

    def set_casedir(self,directory):
        assert os.path.isdir(directory), "not a directory"
        self.casedir = os.path.abspath(directory)

    def set_mesh(self, name, path_to_mesh):
        abspath = os.path.abspath(path_to_mesh)
        self.mesh_dict[name] = abspath

    def load_mesh_dict(self):
        for k, v in self.mesh_dict.items():
            self.mesh_loaded_dict[k] = load_mesh(v)

    def calc_gradients(self, meshname, arr_name):
        assert meshname in self.mesh_loaded_dict.keys(), "can't find mesh in mesh_loaded_dict. run set_mesh and load_mesh_dict"
        mesh = self.mesh_loaded_dict[meshname]
        self.mesh_loaded_dict[meshname] = mesh_scalar_gradients(mesh, self.var_dict[arr_name])

class CascadeCase(AbstractCase):
    def __init__(self, name, vartype, simtype):
        super().__init__(name, vartype, simtype)
        self.x_pos = {"x_pos1": None,
                      "x_pos2": None}
        self.CascadeCoeffs = CascadeCoeffs()

        self.midspan_z = None
        self.machine_type = None


    def set_machine_type(self, mtype):
        allowed = ["compressor", "turbine"]
        assert mtype in allowed, "machine type not allowed. choose between compressor or turbine"
        self.machine_type = mtype

    def set_x_pos(self, x_pos1, x_pos2):

        assert type(x_pos1) == float, "x_pos1 needs to be a float"
        assert type(x_pos2) == float, "x_pos2 needs to be a float"

        bounds = self.mesh_loaded_dict["fluid"].bounds
        x_1 = bounds[0]
        x_2 = bounds[1]

        assert x_pos1 > x_1, "x_pos1 out of fluid-domain-bounds"
        assert x_pos2 < x_2, "x_pos1 out of fluid-domain-bounds"

        self.x_pos["x_pos1"] = x_pos1
        self.x_pos["x_pos2"] = x_pos2

    def get_midspan_z(self):
        if not self.midspan_z:
            mesh = self.mesh_loaded_dict["fluid"]
            self.midspan_z, z = slice_midspan_z(mesh)
        return self.midspan_z, z

case_types = {"CascadeCase" : CascadeCase}
