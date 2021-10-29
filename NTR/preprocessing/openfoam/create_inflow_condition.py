import os
from scipy import interpolate
import numpy as np


from NTR.utils.pyvista_utils import load_mesh


def create_pointfile_txt(points):
    point_header = r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  plus                                  |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant/boundaryData/INLET";
    object      points;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"""

    txt = point_header
    txt += "\n("
    for pt in points:
        txt += "\n("
        txt += str(pt[0])
        txt += "\t\t"
        txt += str(pt[1])
        txt += "\t\t"
        txt += str(pt[2])
        txt += ")"
    txt += "\n)"
    return txt


def create_inflow_condition(simcase_directory, inlet_vtk_mesh_path, onedimensional_data_dicts, positions):
    dirstructure = ["constant", "boundaryData", "INLET", "0"]
    inlet = load_mesh(inlet_vtk_mesh_path)
    inlet.clear_arrays()

    for dir_idx in range(len(dirstructure)):
        dirs = dirstructure[0:dir_idx + 1]
        if not os.path.isdir(os.path.join(simcase_directory, *dirs)):
            os.mkdir(os.path.join(simcase_directory, *dirs))

        if dirs[-1] == "INLET":
            pointfiletxt = create_pointfile_txt(inlet.points)
            pointfilename = "points"
            with open(os.path.join(simcase_directory, *dirs, pointfilename), "w") as fobj:
                fobj.write(pointfiletxt)

        elif dirs[-1] == "0":
            create_boundary_data(inlet.points, onedimensional_data_dicts, os.path.join(simcase_directory, *dirs),
                                 positions)


def create_boundary_data(points, onedimensional_data_dicts, directory, positions):
    for variable_name, values in onedimensional_data_dicts.items():
        datatext = ""
        datatext += "(\n"
        if np.isscalar(values[0]):
            for point in points:

                f = interpolate.interp1d(positions, values)
                ynew = f(point[1])  # use interpolation function returned by `interp1d`
                datatext += str(ynew)
                datatext += "\n"
        else:
            no_columes = values.shape[1]
            for point in points:
                datatext += "( "
                for idx in range(no_columes):
                    column_vals = values[:,idx]
                    f = interpolate.interp1d(positions, column_vals)
                    ynew = f(point[1])  # use interpolation function returned by `interp1d`
                    datatext += str(ynew)
                    datatext += "\t"
                datatext += " )\n"
        datatext += ")\n"
        with open(os.path.join(directory, variable_name), "w") as fobj:
            fobj.write(datatext)
