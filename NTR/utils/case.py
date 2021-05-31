

class abstract_case:
    def __init__(self,name):
        self.name = name


class cascade_case(abstract_case):
    def __init__(self,name):
        super().__init__(name)

        self.mesh_dict = {"fluid":"",
                          "blade":"",
                          "pitch_periodic":"",
                          "span_periodic":"",
                          }
