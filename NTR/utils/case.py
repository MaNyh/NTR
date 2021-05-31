import os


class abstract_case:
    def __init__(self, name):
        self.name = name
        self.mesh_dict = {}

    def set_mesh(self, name, path_to_mesh):
        abspath = os.path.abspath(path_to_mesh)
        self.mesh_dict[name] = abspath


class cascade_case(abstract_case):
    def __init__(self, name):
        super().__init__(name)


