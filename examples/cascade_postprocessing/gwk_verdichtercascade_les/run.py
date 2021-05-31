from NTR.utils.case import cascade_case
from NTR.postprocessing.createProfileData import createProfileData


case = cascade_case("GWK_Verdichterkaskade_LES")

case.set_mesh("fluid", "./VTK/06_GWKVD_81000.vtk")
case.set_mesh("blade", "./VTK/BLADE/BLADE_81000.vtk")

createProfileData(case.mesh_dict["fluid"])
