from NTR.preprocessing.openfoam.create_probes import createProbesProfileDict, createProbesStreamlineDict
from NTR.utils.case import CascadeCase

case = CascadeCase("GWK_Verdichterkaskade_LES")

case.set_casedir(".")

case.set_mesh("fluid", "ressources/VTK/05_GWKVD_109000.vtk")
case.set_mesh("blade", "ressources/VTK/BLADE/BLADE_109000.vtk")
case.load_mesh_dict()

case.CascadeCoeffs.set_alpha(0.01)

createProbesProfileDict(case, 4, 1, 1, ".", tolerance=1e-6)
createProbesStreamlineDict(case, 20, ".", 1, 135, 115, 0.0765)
