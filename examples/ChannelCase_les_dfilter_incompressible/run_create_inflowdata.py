from NTR.preprocessing.openfoam.create_inflow_condition import create_inflow_condition
import numpy as np


def create_onedimensional_datadict():
    positions = np.arange(0, 3, 0.01)
    vars = ["L", "U"]
    datadict = {}
    for var in vars:
        if var == "L":
            datadict[var] = positions
        elif var == "U":
            datadict[var] = np.stack([positions, positions, positions]).T
    return datadict, positions


datadict, positions = create_onedimensional_datadict()
create_inflow_condition(r"D:\CodingProjects\NTR\examples\ChannelCase_les\02_Simcase",
                        r"D:\CodingProjects\NTR\examples\ChannelCase_les\03_Solution\mittelung99Domain_INLET_148000.vtk",
                        datadict, positions)
