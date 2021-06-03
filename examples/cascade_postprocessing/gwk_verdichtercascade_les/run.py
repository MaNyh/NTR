from NTR.utils.case import CascadeCase
from NTR.postprocessing.createProfileData import createProfileData


case = CascadeCase("GWK_Verdichterkaskade_LES")

case.set_mesh("fluid", "./VTK/06_GWKVD_81000.vtk")
case.set_mesh("blade", "./VTK/BLADE/BLADE_81000.vtk")

case.set_x_pos(-0.05, 0.2)

case.fluid_coeffs.set_kappa(1.4)
case.fluid_coeffs.set_M(28.86e-3)
case.fluid_coeffs.set_mu(20e-6)
case.fluid_coeffs.set_p_k(8186.0)  # p_k=23260.0
case.fluid_coeffs.set_kappa(1.4)
case.fluid_coeffs.set_As(1.458e-06)
case.fluid_coeffs.set_cp(1004.5)
case.fluid_coeffs.set_l(0.0927)  # l=0.10
case.fluid_coeffs.set_Ts(110.4)

createProfileData(case)
