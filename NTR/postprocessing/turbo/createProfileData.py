import matplotlib.pyplot as plt
import os
import numpy as np
import math

from NTR.utils.geom_functions.geom_utils import GetProfileValuesMidspan, getPitchValuesB2BSliceComplete
from NTR.utils.fluid_functions.aeroFunctions import Ma, Ma_is, Ma_is_x, Re, Re_is, p_t_is, T_t_is, AVDR, Beta, \
    calcPos2ValuesByAmecke
from NTR.utils.fluid_functions.thermoFunctions import Sutherland_Law
from NTR.utils.mathfunctions import absvec_array
from NTR.utils.simFunctions import sort_value2, sort_value3
from NTR.utils.externals.tecplot.tecplot_functions import writeTecplot1DFile
from NTR.postprocessing.turbo.profile_loading import calc_inflow_cp


def createProfileData(mesh, midspan_z, alpha, post_slice_1_x, post_slice_2_x, output_path, kappa, R_L, p_k, As, l, cp,
                      Ts):
    """

    :param mesh: volume mesh
    :param midspan_z: [m]
    :param alpha: dimensionless parameter for geo-extraction
    :param post_slice_1_x: [m]
    :param post_slice_2_x: [m]
    :param output_path: os.path-like
    :param kappa: isentropic exponent
    :param R_L: specific gas constant
    :param p_k: ???
    :param As:
    :param l: characteristic length
    :param cp: heat capacity isobar
    :return:
    """
    values_ss, values_ps = GetProfileValuesMidspan(mesh, alpha, midspan_z)

    x_ss = values_ss["x_ss"]
    y_ss = values_ss["y_ss"]
    p_ss = values_ss["p"]

    x_ps = values_ps["x_ps"]
    y_ps = values_ps["y_ps"]
    p_ps = values_ps["p"]

    #plt.ioff()
    plt.figure(figsize=(8, 8))
    plt.plot(x_ss, y_ss, '-r', lw=1)
    plt.plot(x_ps, y_ps, '-b', lw=1)

    plt.axis('equal')
    plt.grid()
    plt.savefig(os.path.join(output_path, 'kontrollplot_profil.pdf'))

    inte_mag_u1, inte_ux1, inte_uy1, inte_uz1, inte_rho1, inte_T1, inte_p1, inte_p_tot1, inte_T_tot1 = calcPostSliceValues(
        mesh, output_path, post_slice_1_x, 1, kappa, R_L)
    inte_mag_u2, inte_ux2, inte_uy2, inte_uz2, inte_rho2, inte_T2, inte_p2, inte_p_tot2, inte_T_tot2 = calcPostSliceValues(
        mesh, output_path, post_slice_2_x, 2, kappa, R_L)

    zeta = (inte_p_tot1 - inte_p_tot2) / (inte_p_tot1 - p_k)

    Ma1 = Ma(inte_mag_u1, kappa, R_L, inte_T1)
    Ma2 = Ma(inte_mag_u2, kappa, R_L, inte_T2)

    Ma_is_2 = Ma_is(p_k, kappa, inte_p1, inte_rho1, inte_mag_u1, R_L, inte_T_tot1)

    Re_is_2 = Re_is(kappa, R_L, l, As, Ma_is_2, p_k, inte_T_tot1, inte_mag_u1, cp, Ts)

    beta1 = 90.0 + math.atan(inte_uy1 / inte_ux1) / 2.0 / math.pi * 360.0
    beta2 = 90.0 + math.atan(inte_uy2 / inte_ux2) / 2.0 / math.pi * 360.0

    nu1 = Sutherland_Law(inte_T1, As, Ts)
    nu2 = Sutherland_Law(inte_T2, As, Ts)

    Re1 = Re(inte_rho1, inte_mag_u1, l, nu1)
    Re2 = Re(inte_rho2, inte_mag_u2, l, nu2)

    AVDR_value = AVDR(inte_rho1, inte_mag_u1, beta1, inte_rho2, inte_mag_u2, beta2)
    delta_beta = beta1 - beta2
    delta_p_static = (inte_p1 - inte_p2) / (inte_p_tot1 - p_k)

    y, array_names, values = getPitchValuesB2BSliceComplete(mesh, post_slice_2_x)

    # brechnung der amecke werte
    p_2_y = values[array_names.index('p')]
    pt_2_y = []

    Ux_2_y = np.asarray(values[array_names.index('U')])[::, 0]
    Uy_2_y = np.asarray(values[array_names.index('U')])[::, 1]

    Mag_U_2_y = absvec_array(values[array_names.index("U")])  # values[array_names.index('Mag_U')]

    T_2_y = values[array_names.index('T')]
    beta_2_y = []
    for i in range(len(p_2_y)):
        beta_2_y.append(Beta(Ux_2_y[i], Uy_2_y[i]))
        pt_2_y.append(p_t_is(kappa, Ma(Mag_U_2_y[i], kappa, R_L, T_2_y[i]), p_2_y[i]))

    Ma2_amecke, beta2_amecke, pt2_amecke, p2_amecke = calcPos2ValuesByAmecke(pt_2_y, beta_2_y, p_2_y, y, inte_p_tot1,
                                                                             kappa=kappa)

    zeta_amecke = (inte_p_tot1 - pt2_amecke) / (inte_p_tot1 - p_k)
    ma_is_amecke = Ma_is_x(kappa, p2_amecke, pt2_amecke)

    x_ss, y_ss, x_zu_l_ax_ss, p_ss, cp_ss, cp_max_ss, ma_is_x_ss, x_ps, y_ps, x_zu_l_ax_ps, p_ps, cp_ps, cp_max_ps, ma_is_x_ps = calcProfileValues(
        p_ss, p_ps, x_ss, inte_p_tot1, output_path, x_ps, y_ss, y_ps, inte_mag_u1, inte_rho1, inte_p1, kappa)


    ans = {}
    ans["zeta"]=zeta
    ans["inte_p1"]=inte_p1
    ans["inte_p_tot1"]=inte_p_tot1
    ans["inte_mag_u1"]=inte_mag_u1
    ans["inte_ux1"]=inte_ux1
    ans["inte_uy1"]=inte_uy1
    ans["inte_uz1"]=inte_uz1
    ans["inte_rho1"]=inte_rho1
    ans["inte_T1"]=inte_T1
    ans["inte_T_tot1"]=inte_T_tot1
    ans["beta1"]=beta1
    ans["nu1"]=nu1
    ans["Re1"]=Re1
    ans["Ma1"]=Ma1
    ans["Ma2"]=Ma2
    ans["inte_p2"]=inte_p2
    ans["inte_p_tot2"]=inte_p_tot2
    ans["inte_mag_u2"]=inte_mag_u2
    ans["inte_ux2"]=inte_ux2
    ans["inte_uy2"]=inte_uy2
    ans["inte_uz2"]=inte_uz2
    ans["inte_rho2"]=inte_rho2
    ans["inte_T2"]=inte_T2
    ans["inte_T_tot2"]=inte_T_tot2
    ans["beta2"]=beta2
    ans["nu2"]=nu2
    ans["Re2"]=Re2
    ans["Ma_is_2"]=Ma_is_2
    ans["Re_is_2"]=Re_is_2
    ans["AVDR_value"]=AVDR_value
    ans["delta_beta"]=delta_beta
    ans["delta_p_static"]=delta_p_static
    ans["Ma2_amecke"]=Ma2_amecke
    ans["beta2_amecke"]=beta2_amecke
    ans["pt2_amecke"]=pt2_amecke
    ans["p2_amecke"]=p2_amecke
    ans["zeta_amecke"]=zeta_amecke
    ans["ma_is_amecke"]=ma_is_amecke
    ans["profileData"] = {"x_zu_l_ax_ss":x_zu_l_ax_ss, "cp_ss":cp_ss, "x_zu_l_ax_ps":x_zu_l_ax_ps, "cp_ps":cp_ps}


    return ans#x_zu_l_ax_ss, cp_ss, x_zu_l_ax_ps, cp_ps


def calcProfileValues(p_ss, p_ps, x_ss, inte_p_tot1, output_path, x_ps, y_ss, y_ps, inte_mag_u1, inte_rho1, inte_p1, kappa):

    p = p_ss + p_ps[::-1]
    p_max = max(p)
    p_te = p_ps[-1]

    cp_ss = []
    cp_max_ss = []
    ma_is_x_ss = []

    x_zu_l_ax_ss = []

    for i in range(len(x_ss)):
        cp_ss.append(calc_inflow_cp(p_ss[i],inte_p_tot1,inte_p1))
        cp_max_ss.append((p_ss[i] - p_te) / (p_max - p_te))
        x_zu_l_ax_ss.append((x_ss[i] - min(x_ss)) / (max(x_ss) - min(x_ss)))
        ma_is_x_ss.append(Ma_is_x(kappa, p_ss[i], inte_p_tot1))

    cp_ps = []
    cp_max_ps = []
    ma_is_x_ps = []

    x_zu_l_ax_ps = []

    for i in range(len(x_ps)):
        cp_ps.append(calc_inflow_cp(p_ps[i],inte_p_tot1,inte_p1))
        cp_max_ps.append((p_ps[i] - p_te) / (p_max - p_te))
        x_zu_l_ax_ps.append((x_ps[i] - min(x_ps)) / (max(x_ps) - min(x_ps)))
        ma_is_x_ps.append(Ma_is_x(kappa, p_ps[i], inte_p_tot1))

    # ACHTUNG, ps wird gedreht, weil calcConcaveHull im mathematischem Drehsinn die
    cp = cp_ss + cp_ps[::-1]
    cp_max = cp_max_ss + cp_max_ps[::-1]
    x_zu_l_ax = x_zu_l_ax_ss + x_zu_l_ax_ps[::-1]
    ma_is_x = ma_is_x_ss + ma_is_x_ps[::-1]

    plt.figure(figsize=(9, 6))
    plt.plot(x_zu_l_ax, cp, '-r', lw=1, label='cp')
    plt.plot(x_zu_l_ax, cp_max, '-b', lw=1, label='cp_max')
    plt.legend()
    ax1 = plt.gca()
    ax1.set_xlim([0, 1])
    ax1.set_ylim([-1, 1])
    plt.grid()
    plt.savefig(os.path.join(output_path, 'kontrollplot_cp.pdf'))
    plt.close('all')

    plt.figure(figsize=(9, 6))
    plt.plot(x_zu_l_ax, ma_is_x, '-r', lw=1, label='Ma_is_x')
    plt.legend()
    ax1 = plt.gca()
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0, 1])
    plt.grid()
    plt.savefig(os.path.join(output_path, 'kontrollplot_ma_is_x.pdf'))
    plt.close('all')

    return x_ss, y_ss, x_zu_l_ax_ss, p_ss, cp_ss, cp_max_ss, ma_is_x_ss, x_ps, y_ps, x_zu_l_ax_ps, p_ps, cp_ps, cp_max_ps, ma_is_x_ps


def calcPostSliceValues(mesh, output_path, x, ind, kappa, R_L):
    cut_plane = mesh.slice(normal="x", origin=(x, 0, 0))
    points = cut_plane.points
    npts = cut_plane.number_of_points

    xx = np.zeros(npts)
    y = np.zeros(npts)
    zz = np.zeros(npts)

    for ii in range(npts):
        pt = points[ii]
        y[ii] = pt[1]
        xx[ii] = pt[0]
        zz[ii] = pt[2]

    mag_u_array = absvec_array(cut_plane.point_arrays["U"])

    nvls = len(mag_u_array)

    ux_array = cut_plane.point_arrays["U"][::, 0]
    uy_array = cut_plane.point_arrays["U"][::, 1]
    uz_array = cut_plane.point_arrays["U"][::, 2]

    rho_array = cut_plane.point_arrays["rho"]
    T_array = cut_plane.point_arrays["T"]
    p_array = cut_plane.point_arrays["p"]

    mag_u = []

    ux = []
    uy = []
    uz = []

    rho = []
    T = []
    p = []
    ma = []

    for i in range(nvls):
        mag_u.append(mag_u_array[i])
        ux.append(ux_array[i])
        uy.append(uy_array[i])
        uz.append(uz_array[i])
        rho.append(rho_array[i])
        T.append(T_array[i])
        p.append(p_array[i])
        ma.append(Ma(mag_u[-1], kappa, R_L, T[-1]))

    xx = sort_value3(y, xx)
    zz = sort_value3(y, zz)
    mag_u = sort_value3(y, mag_u)

    ux = sort_value3(y, ux)
    uy = sort_value3(y, uy)
    uz = sort_value3(y, uz)

    rho = sort_value3(y, rho)
    T = sort_value3(y, T)
    p = sort_value3(y, p)
    y, ma = sort_value2(y, ma)
    # totaldruck berechnen
    p_tot = []
    T_tot = []
    for i in range(len(p)):
        # kappa ma
        p_tot.append(p_t_is(kappa, ma[i], p[i]))
        T_tot.append(T_t_is(kappa, ma[i], T[i]))

    inte_ux = mass_average(y, ux, rho, ux)
    inte_uy = mass_average(y, uy, rho, ux)
    inte_uz = mass_average(y, uz, rho, ux)
    inte_mag_u = np.sqrt(inte_ux ** 2 + inte_uy ** 2 + inte_uz ** 2)

    inte_rho = area_average(y, rho)
    inte_T = area_average(y, T)
    inte_p = area_average(y, p)

    inte_p_tot = mass_average(y, p_tot, rho, ux)
    inte_T_tot = mass_average(y, T_tot, rho, ux)

    f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(1, 7, sharey=True)

    ax1.plot(mag_u, y)
    # ax1.plot(mean_mag_u,inte_mag_u)
    ax1.set_title('Mag U')

    ax2.plot(ux, y)
    ax2.set_title('Ux')

    ax3.plot(uy, y)
    ax3.set_title('Uy')

    ax4.plot(uz, y)
    ax4.set_title('Uz')

    ax5.plot(rho, y)
    ax5.set_title('rho')

    ax6.plot(T, y)
    ax6.set_title('T')

    ax7.plot(p, y)
    ax7.set_title('p')

    plt.grid()
    # plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'kontrollplot_auswerteebene_' + str(ind) + '.pdf'))

    values = [[xx, y, zz, mag_u, ux, uy, uz, p, rho, T, ma, T_tot, p_tot]]

    data_file_path = os.path.join(output_path, 'values_auswerteebene_' + str(ind) + '.dat')
    writeTecplot1DFile(data_file_path, ['X', 'Y', 'Z', 'mag_U', 'U', 'V', 'W', 'p', 'rho', 'T', 'Ma', 'Tt', 'pt'],
                       ['values_auswerteebene_' + str(ind)], values, 'Profildaten')

    return inte_mag_u, inte_ux, inte_uy, inte_uz, inte_rho, inte_T, inte_p, inte_p_tot, inte_T_tot


def area_average(y, var):
    area_average_val = np.trapz(var, x=y) / (max(y) - min(y))

    return area_average_val


def mass_average(y, var, rho, velo):
    mass_flow_local = []
    mass_flow_dot_var = []

    for i in range(len(y)):
        mass_flow_local.append(rho[i] * velo[i])
        mass_flow_dot_var.append(rho[i] * velo[i] * var[i])

    mass_flow_global = np.trapz(mass_flow_local, x=y)

    mass_flow_dot_var_int = np.trapz(mass_flow_dot_var, x=y)

    mass_average_val = mass_flow_dot_var_int / mass_flow_global

    return mass_average_val
