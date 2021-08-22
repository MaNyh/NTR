# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 00:27:21 2018

@author: Mark
"""
import numpy as np


def calcCf(tau_w, U_inf, rho_inf):
    # berechnet den reibungsbeiwert

    # tau_w: Wandschubspannung
    # U_inf: Zuestroemgeschwindigkeit
    # rho_inf: Zuestroemdichte

    cf = 2 * tau_w / (rho_inf * U_inf ** 2)
    return cf



def calcWallShearStress(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, face_normal, rho, nu, p):
    # face_normal: normal des faces aufm blade, zeigt vom blade weg
    # nu: kinematische viskositaet

    grad_U_tensor = np.array([[dudx, dudy, dudz], [dvdx, dvdy, dvdz], [dwdx, dwdy, dwdz]])
    I = np.identity(3)
    mu = rho * nu
    mu = nu
    mean_tau = -p * I + mu * (
            grad_U_tensor + grad_U_tensor.transpose() - (2 / 3.0) * grad_U_tensor.transpose().trace() * I)
    # print(mean_tau)
    shear_stress_at_wall = np.dot(mean_tau, face_normal)
    wall_normal_shear_stress = np.dot(shear_stress_at_wall, face_normal)
    wall_tangent_shear_stress = shear_stress_at_wall - wall_normal_shear_stress * face_normal
    return wall_tangent_shear_stress

