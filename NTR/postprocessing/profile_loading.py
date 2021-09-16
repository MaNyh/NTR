import os

from NTR.utils.mathfunctions import vecAbs
from NTR.database.case_structure import case_dirs

class profile_loading:
    def __init__(pathsettings_yml):
        self.ssPoly = None
        self.psPoly = None

        self.case_path = os.path.abspath(os.path.dirname(settings_yaml))
        self.settings = yaml_dict_read(settings_yaml)

        self.solutionpath = os.path.join(case_path, case_dirs["solution"])
        self.datpath = os.path.join(case_path, case_dirs["data"])

        assert os.path.isdir(self.datpath), "no data-directory found. either nothing was created or it got deleted."
        assert os.path.isdir(self.solutionpath), "no solution-directory found."

    def calc_loading_from_volmesh(self):
        meshname = self.settings[]

    def calc_x_c(self):
        x_ss = self.ssPoly.points[::,0]
        x_ps = self.psPoly.points[::,0]
        c_ss = vecAbs(self.ssPoly.points[-1]-self.ssPoly.points[0])
        c_ps = vecAbs(self.psPoly.points[-1]-self.psPoly.points[0])
        assert c_ss==c_ps, "camberlength from side-polys not equal"
        self.ssPoly["x_c"]=x_ss/c_ss
        self.psPoly["x_c"]=x_ps/c_ps


def calc_cp(p_blade, p_inlet, c_inlet, rho):
    return p_blade-p_inlet/(1/2*rho*c_inlet**2)
