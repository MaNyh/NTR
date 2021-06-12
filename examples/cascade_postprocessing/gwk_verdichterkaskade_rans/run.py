from NTR.utils.case import CascadeCase
from NTR.postprocessing.createProfileData import createProfileData


case = CascadeCase("GWK_Verdichterkaskade_RANS", "openfoam_instantanious")
case.set_machine_type("compressor")
case.set_casedir(".")

case.set_mesh("fluid", "./VTK/06_GWKVD_RANS_22900.vtk")
case.set_mesh("blade", "./VTK/BLADE/BLADE_22900.vtk")
case.load_mesh_dict()

case.calc_gradients("fluid", "U")

case.set_x_pos(-0.08, 0.2)

case.CascadeCoeffs.set_alpha(0.01)

case.FluidCoeffs.set_kappa(1.4)
case.FluidCoeffs.set_M(28.86e-3)
case.FluidCoeffs.set_mu(20e-6)
case.FluidCoeffs.set_As(1.458e-06)
case.FluidCoeffs.set_cp(1004.5)
case.FluidCoeffs.set_Ts(110.4)

case.FluidCoeffs.set_l(0.13)  # l=0.10
case.FluidCoeffs.set_p_k(95000.0)  # Austrittsdruck

createProfileData(case)
