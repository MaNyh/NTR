from NTR.preprocessing.openfoam.create_probes import createProbesProfileDict, createProbesStreamlineDict
from NTR.utils.case import CascadeCase

case = CascadeCase("GWK_Verdichterkaskade_LES","openfoam_timeaveraged", "placeholder_simtype")

case.set_casedir(".")

case.set_mesh("fluid", "vtk_results/VTK/06_GWKVD_317000.vtk")
case.set_mesh("blade", "vtk_results/VTK/BLADE/BLADE_317000.vtk")
case.load_mesh_dict()

case.CascadeCoeffs.set_alpha(0.01)

createProbesProfileDict(case, 20, 20, 1, ".", tolerance=1e-6)
createProbesStreamlineDict(case, 20, ".", 1, 135, 115, 0.0765)
