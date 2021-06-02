from NTR.utils.case import cascade_case
from NTR.postprocessing.createProfileData import createProfileData


case = cascade_case("GWK_Verdichterkaskade_LES")

case.set_mesh("fluid", "./VTK/06_GWKVD_81000.vtk")
case.set_mesh("blade", "./VTK/BLADE/BLADE_81000.vtk")

case.set_x_pos(-0.05, 0.25)

case.fluid_coeffs.set_kappa(1.4)
case.fluid_coeffs.set_M(28.86e-3)
case.fluid_coeffs.set_mu(20e-6)

createProfileData(case)
