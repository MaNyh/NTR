from NTR.utils.case import CascadeCase
from NTR.postprocessing.createProfileData import createProfileData


case = CascadeCase("GWK_Verdichterkaskade_LES")
case.set_casedir(".")

case.set_mesh("fluid", "./VTK/06_GWKVD_81000.vtk")
#case.set_mesh("blade", "./VTK/BLADE/BLADE_81000.vtk")
case.load_mesh_dict()

case.calc_gradients("fluid", "U")

case.set_x_pos(-0.05, 0.2)

case.CascadeCoeffs.set_alpha(0.01)

case.FluidCoeffs.set_kappa(1.4)
case.FluidCoeffs.set_M(28.86e-3)
case.FluidCoeffs.set_mu(20e-6)
case.FluidCoeffs.set_p_k(23260.0)  # p_k=23260.0 8186.0
case.FluidCoeffs.set_kappa(1.4)
case.FluidCoeffs.set_As(1.458e-06)
case.FluidCoeffs.set_cp(1004.5)
case.FluidCoeffs.set_l(0.0927)  # l=0.10
case.FluidCoeffs.set_Ts(110.4)

createProfileData(case)
