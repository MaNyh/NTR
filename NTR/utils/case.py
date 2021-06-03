import os

from NTR.utils.coefficients import FluidCoeffs


class AbstractCase:
    def __init__(self, name):
        self.name = name
        self.mesh_dict = {}
        self.fluid_coeffs = FluidCoeffs()

    def set_mesh(self, name, path_to_mesh):
        abspath = os.path.abspath(path_to_mesh)
        self.mesh_dict[name] = abspath


class CascadeCase(AbstractCase):
    def __init__(self, name):
        super().__init__(name)
        self.x_pos = {"x_pos1": None,
                      "x_pos2": None}

    def set_x_pos(self, x_pos1, x_pos2):

        assert type(x_pos1) == float, "x_pos1 needs to be a float"
        assert type(x_pos2) == float, "x_pos2 needs to be a float"

        self.x_pos["x_pos1"] = x_pos1
        self.x_pos["x_pos2"] = x_pos2
