# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 20:46:24 2018

@author: ziesse
"""
import numpy as np


# Totaldruckverlustbeiwert
def calcZeta(pt1, pt2x, p2):
    zeta = (pt1 - pt2x) / (pt1 - p2)
    return zeta


# lokale isentrope Mach-Zahl
def Ma_is_x(kappa, p, pt):
    if p > 0 and pt > 0 and pt >= p:
        y = np.sqrt(2 / (kappa - 1) * ((pt / p) ** ((kappa - 1) / kappa) - 1))
    else:
        y = 0.0
    return y


# Druckbeiwert
def calcCp(px, pt1, p2):
    cp = (px - p2) / (pt1 - p2)
    return cp


# isentrope reynoldszahl
def Re_is(k, R, l_chord, beta_s, Ma2th, pk, T1, Mag_U, cp, S):
    Tt1 = T_t_is(k, Ma(Mag_U, k, R, T1), T1)
    t_iso = Tt1 / (1.0 + ((k - 1.0) / 2.0) * pow(Ma2th, 2))
    y = np.sqrt(k / R) * l_chord / beta_s * (Ma2th * pk * (t_iso + S)) / pow(t_iso, 2)
    return y


# isenstrope machzahl
def Ma_is(pk, kappa, p1, rho1, Mag_U, R, T1):
    pt1 = p_t_is(kappa, Ma(Mag_U, kappa, R, T1), p1)
    q2th = pt1 - pk
    y = np.sqrt(2.0 / (kappa - 1.0) * (pow(1.0 + (q2th / pk), (kappa - 1.0) / kappa) - 1.0))
    return y


def calcBPValues(p_in, rho_in, T_in, Ux_in, Uy_in, p_out, l_chord, kappa=1.4, R=287.058, beta_s=1.458e-06, S=110.4):
    Mag_U = np.sqrt(Ux_in ** 2 + Uy_in ** 2)
    pt_in = p_t_is(kappa, Ma(Mag_U, kappa, R, T_in), p_in)
    print(pt_in)
    q2th = pt_in - p_out
    Ma_is_2 = np.sqrt(2.0 / (kappa - 1.0) * (pow(1.0 + (q2th / p_out), (kappa - 1.0) / kappa) - 1.0))
    print(Ma_is_2)

    Tt_in = T_t_is(kappa, Ma(Mag_U, kappa, R, T_in), T_in)
    t_iso = Tt_in / (1.0 + ((kappa - 1.0) / 2.0) * pow(Ma_is_2, 2))
    Re_is_2 = np.sqrt(kappa / R) * l_chord / beta_s * (Ma_is_2 * p_out * (t_iso + S)) / pow(t_iso, 2)
    print(Re_is_2)


def Ma(c, kappa, R_L, T):
    Ma = c / ((kappa * R_L * T) ** (0.5))
    return Ma


def AVDR(rho_1, mag_u_1, beta_1, rho_2, mag_u_2, beta_2):
    AVDR = rho_2 * mag_u_2 * np.sin(np.deg2rad(beta_2)) / (rho_1 * mag_u_1 * np.sin(np.deg2rad(beta_1)))
    return AVDR


def p_t_is(kappa, ma, p):
    # https://www.grc.nasa.gov/www/BGH/isentrop.html
    p_t_is = p * pow(1.0 + (kappa - 1.0) / 2.0 * pow(ma, 2.0), (kappa / (kappa - 1.0)))

    return p_t_is


def p_is(kappa, ma, p_t_is):
    # https://www.grc.nasa.gov/www/BGH/isentrop.html
    p_is = p_t_is / pow(1.0 + (kappa - 1.0) / 2.0 * pow(ma, 2.0), (kappa / (kappa - 1.0)))

    return p_is


def T_t_is(kappa, ma, T):
    # https://www.grc.nasa.gov/www/BGH/isentrop.html
    T_t_is = T / (((1.0 + (kappa - 1.0) * 0.5 * ma ** 2.0)) ** (-1.0))

    return T_t_is


def T_is(kappa, ma, Tt):
    # https://www.grc.nasa.gov/www/BGH/isentrop.html
    T = Tt / (1 + ((kappa - 1) / 2.0) * ma ** 2)

    return T


def Re(rho, velo, d, nu):
    Re = rho * velo * d / (nu)
    return Re


def Mass_Average(y, var, rho, velo):
    # nachfolgende Funktion massenstrom mittelt eine 1d-linie

    mass_flow_local = []
    mass_flow_dot_var = []

    for i in range(len(y)):
        mass_flow_local.append(rho[i] * velo[i])
        mass_flow_dot_var.append(rho[i] * velo[i] * var[i])

    mass_flow_global = np.trapz(mass_flow_local, x=y)

    mass_flow_dot_var_int = np.trapz(mass_flow_dot_var, x=y)

    mass_average_val = mass_flow_dot_var_int / mass_flow_global

    return mass_average_val


def Flux_Average(y, var, velo):
    # nachfolgende Funktion flussmittelt eine 1d-linie

    flux_local = []
    flux_var = []

    for i in range(len(y)):
        flux_local.append(velo[i])
        flux_var.append(velo[i] * var[i])

    flux_flow_global = np.trapz(flux_local, x=y)

    flux_flow_dot_var_int = np.trapz(flux_var, x=y)

    flux_average_val = flux_flow_dot_var_int / flux_flow_global

    return flux_average_val


def Beta(Ux, Uy):
    beta = 90.0 + np.rad2deg(np.arctan(Uy / float(Ux)))

    return beta


def Area_Average(y, var):
    # nachfolgende Funktion flaechenmittelt eine 1d-linie

    area_average_val = np.trapz(var, x=y) / (max(y) - min(y))

    return area_average_val


def Sr(U_wake, t_wake, l_ax, c_ax):
    # Berechnet die Strouhal-Zahl
    # U-Wake:   Rotations-/Translationsgeschwindigkeit der Wakes
    # t_wake:   Teilung der Wakes
    # l_ax:     axiale Sehnenlänge des zu untersuchenden Profils
    # c_ax:     axiale Zustroemgeschwindigkeit des Profils
    Sr = (l_ax / t_wake) * (U_wake / c_ax)
    return Sr


def Sr2(U_wake, t_wake, l, U_abs_2):
    # Berechnet die Strouhal-Zahl
    # U-Wake:   Rotations-/Translationsgeschwindigkeit der Wakes
    # t_wake:   Teilung der Wakes
    # l:        Sehnenlaenge
    # U_abs_2:     Absolute Abstroemgeschwindigkeit
    Sr = (U_wake / t_wake) * (l / U_abs_2)
    return Sr


def calcPos2ValuesByAmecke(pt_2_y, beta_2_y, p_2_y, y, pt_1, kappa=1.4):
    # winkel umerechnen in amecke format

    def calcQ(kappa, p_zu_pt):
        q_zu_pt = (kappa / (kappa - 1.0)) * (p_zu_pt) ** (1.0 / kappa) * \
                  (1 - (p_zu_pt) ** ((kappa - 1.0) / kappa))

        return q_zu_pt

    def calcTheta(kappa, p_zu_pt):
        theta = np.sqrt(2.0 / (kappa - 1.0)) * np.sqrt(((kappa + 1.0) / 2.0) ** ((kappa + 1.0) / (kappa - 1.0))) * \
                p_zu_pt ** (1 / kappa) * \
                np.sqrt(1 - p_zu_pt ** ((kappa - 1.0) / kappa))

        return theta

    y_rel = []

    for i in range(len(y)):
        y_rel.append((y[i] - min(y)) / (max(y) - min(y)))

    def calcI1(pt_2_y, beta_2_y, p_2_y, y_rel, pt_1, kappa):

        values = []

        for i in range(len(pt_2_y)):
            theta = calcTheta(kappa, p_2_y[i] / pt_2_y[i])

            values.append((pt_2_y[i] / pt_1) * theta * np.sin(np.deg2rad(beta_2_y[i])))

        I1 = np.trapz(values, x=y_rel)
        return I1

    def calcI2(pt_2_y, beta_2_y, p_2_y, y_rel, pt_1, kappa):

        values = []

        for i in range(len(pt_2_y)):
            q_zu_p = calcQ(kappa, p_2_y[i] / pt_2_y[i])
            values.append((pt_2_y[i] / pt_1) *
                          (2 * q_zu_p * np.sin(np.deg2rad(beta_2_y[i])) * np.sin(np.deg2rad(beta_2_y[i])) + p_2_y[i] /
                           pt_2_y[i]))

        I2 = np.trapz(values, x=y_rel)

        return I2

    def calcI3(pt_2_y, beta_2_y, p_2_y, y_rel, pt_1, kappa):

        values = []

        for i in range(len(pt_2_y)):
            q_zu_p = calcQ(kappa, p_2_y[i] / pt_2_y[i])
            values.append(2 * (pt_2_y[i] / pt_1) * q_zu_p *
                          np.sin(np.deg2rad(beta_2_y[i])) *
                          np.cos(np.deg2rad(beta_2_y[i])))

        I3 = np.trapz(values, x=y_rel)

        return I3

    def calcMa2(kappa, I1, I2, I3):

        ma = np.sqrt(((kappa + 1.0) / 2.0) ** (2.0 / (kappa - 1.0)) * (I2 ** 2.0 / I1 ** 2) * (
                0.5 - (2.0 / (kappa + 1.0)) ** (2.0 / (kappa - 1)) * (I1 ** 2 / I2 ** 2) + (
                    (kappa + 1.0) / (2.0 * kappa)) * (I3 ** 2 / I2 ** 2) -
                np.sqrt(0.25 - (2.0 / (kappa + 1.0)) ** (2.0 / (kappa - 1)) * (I1 ** 2 / I2 ** 2) + (
                        (kappa ** 2 - 1.0) / (4.0 * kappa ** 2)) * (I3 ** 2 / I2 ** 2))))

        return ma

    I1 = calcI1(pt_2_y, beta_2_y, p_2_y, y_rel, pt_1, kappa)
    I2 = calcI2(pt_2_y, beta_2_y, p_2_y, y_rel, pt_1, kappa)
    I3 = calcI3(pt_2_y, beta_2_y, p_2_y, y_rel, pt_1, kappa)
    Ma2 = calcMa2(kappa, I1, I2, I3)

    p2_zu_pt2 = (1.0 - ((kappa - 1.0) / (kappa + 1)) * Ma2 ** 2) ** (kappa / (kappa - 1.0))
    theta2 = calcTheta(kappa, p2_zu_pt2)
    q2_zu_pt2 = calcQ(kappa, p2_zu_pt2)

    beta2 = np.rad2deg(np.arccos((I3 / I1) * (theta2 / (2.0 * q2_zu_pt2))))
    pt2 = pt_1 * (I1 / (theta2 * np.sin(np.deg2rad(beta2))))
    p2 = (1 - ((kappa - 1.0) / (kappa + 1.0)) * Ma2 ** 2.0) ** (kappa / (kappa - 1.0)) * pt2
    beta2 = beta2

    return Ma2, beta2, pt2, p2


def estimateBP(Re_2th_soll, Ma_2th_soll, l_chord, Tt1, pt1, pk, R=287.058, beta_s=1.458e-6, S=110.4, kappa=1.4):
    err = 1.0

    while err > 1e-12:

        q2th = pt1 - pk

        err_Ma = 0.0
        err_Re = 0.0

        try:

            Ma2th = np.sqrt(2.0 / (kappa - 1.0) * (pow(1.0 + (q2th / pk), (kappa - 1.0) / kappa) - 1.0))
        except:
            err_Ma = 100.0

        t_iso = Tt1 / (1.0 + ((kappa - 1.0) / 2.0) * pow(Ma2th, 2))

        try:
            Re2th = np.sqrt(kappa / R) * l_chord / beta_s * (Ma2th * pk * (t_iso + S)) / pow(t_iso, 2)
        except:
            err_Re = 100.0

        print("pt1\t= %s\tpk\t= %s\tMa_2,th\t= %s\tRe_2,th\t= %s" % (pt1, pk, Ma2th, Re2th))

        err_Re += abs(Re2th - Re_2th_soll) / Re_2th_soll

        err_Ma += abs(Ma2th - Ma_2th_soll) / Ma_2th_soll

        err = err_Ma + err_Re

        # sleep(0.01)

        print("err = %s" % (err))

        # if err_Ma<err_Re:
        pk *= 1.0 + max(-0.2, min(0.2, 0.5 * (Re_2th_soll - Re2th) / Re_2th_soll))
        # else:
        q2th *= 1.0 + max(-0.2, min(0.2, 0.5 * (Ma_2th_soll - Ma2th) / Ma_2th_soll))

        pt1 = pk + q2th



