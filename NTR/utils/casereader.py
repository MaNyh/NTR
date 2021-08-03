import os

from NTR.utils.case import case_types
from NTR.utils.filehandling import yaml_dict_read


def case_from_dict(case_yml):
    case_dict = yaml_dict_read(case_yml)
    case_dir = os.path.abspath(os.path.dirname(case_yml))
    sets = case_dict.keys()

    if "case_settings" in sets:
        name = case_dict["case_settings"]["name"]
        ctype = case_dict["case_settings"]["case_type"]
        stype = case_dict["case_settings"]["sim_type"]
        vtype = case_dict["case_settings"]["var_type"]
        mtype = case_dict["case_settings"]["machine_type"]
        case = case_types[ctype](name, vtype, stype)
        case.set_casedir(case_dir)
        case.set_machine_type(mtype)

    if "mesh_dict" in sets:
        meshdict = case_dict["mesh_dict"]
        for mesh_name, mesh_path in meshdict.items():
            mesh_path = os.path.join(case.casedir, mesh_path)
            case.set_mesh(mesh_name, mesh_path)
        case.load_mesh_dict()

    if "post_process_settings" in sets:
        if "calc_gradients" in case_dict["post_process_settings"]:
            mesh_name = case_dict["post_process_settings"]["calc_gradients"]["mesh"]
            scalarfield = case_dict["post_process_settings"]["calc_gradients"]["scalar"]
            case.calc_gradients(mesh_name, scalarfield)

        if "set_xs" in case_dict["post_process_settings"]:
            x1pos = case_dict["post_process_settings"]["set_xs"]["set_x1"]
            x2pos = case_dict["post_process_settings"]["set_xs"]["set_x2"]
            case.set_x_pos(x1pos, x2pos)

        if "alpha" in case_dict["post_process_settings"]:
            alpha = case_dict["post_process_settings"]["alpha"]
            case.CascadeCoeffs.set_alpha(alpha)

    if "fluid_coeffs" in case_dict:
        kappa = case_dict["fluid_coeffs"]["kappa"]
        case.FluidCoeffs.set_kappa(kappa)
        M = case_dict["fluid_coeffs"]["M"]
        case.FluidCoeffs.set_M(M)
        mu = case_dict["fluid_coeffs"]["mu"]
        case.FluidCoeffs.set_mu(mu)
        As = case_dict["fluid_coeffs"]["As"]
        case.FluidCoeffs.set_As(As)
        cp = case_dict["fluid_coeffs"]["cp"]
        case.FluidCoeffs.set_cp(cp)
        Ts = case_dict["fluid_coeffs"]["Ts"]
        case.FluidCoeffs.set_Ts(Ts)
        l = case_dict["fluid_coeffs"]["l"]
        case.FluidCoeffs.set_l(l)
        p_k = case_dict["fluid_coeffs"]["p_k"]
        case.FluidCoeffs.set_p_k(p_k)

    return case
