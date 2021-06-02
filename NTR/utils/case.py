import os

class fluid_coeffs:
    #https://www.chemie.de/lexikon/Universelle_Gaskonstante.html
    R = 8.314472

    def __init__(self):
        self.kappa = None
        self.R_L = None
        self.mu = None
        self.M = None

    def set_kappa(self, kappa):
        assert type(kappa) == float, "kappa needs to be a float"
        self.kappa = kappa

    def set_R_L(self, R_L):
        assert type(R_L) == float, "R_L needs to be a float"
        self.R_L = R_L

    def set_mu(self, mu):
        assert type(mu) == float, "mu needs to be a float"
        self.mu = mu

    def set_M(self, M):
        assert type(M) == float, "M needs to be a float"
        self.M = M
        self.R_L = self.R / self.M

class abstract_case:
    def __init__(self, name):
        self.name = name
        self.mesh_dict = {}
        self.fluid_coeffs = fluid_coeffs()

    def set_mesh(self, name, path_to_mesh):
        abspath = os.path.abspath(path_to_mesh)
        self.mesh_dict[name] = abspath


class cascade_case(abstract_case):
    def __init__(self, name):
        super().__init__(name)
        self.x_pos = {"x_pos1":None,
                      "x_pos2":None}

    def set_x_pos(self,x_pos1,x_pos2):

        assert type(x_pos1) == float, "x_pos1 needs to be a float"
        assert type(x_pos2) == float, "x_pos2 needs to be a float"

        self.x_pos["x_pos1"] = x_pos1
        self.x_pos["x_pos2"] = x_pos2


