#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 07:43:06 2018

@author: ziesse
"""

import os
import vtk

"""
from NTR.Libs.geomFunctions import calc_vk_hk
from NTR.Libs.geomFunctions import getGeom2DVTUSLice
from NTR.Libs.geomFunctions import sortProfilPoints
from NTR.Libs.geomFunctions import rotateTensorVector
from NTR.Libs.boundaryLayerFunctions import calcWallShearStress
from NTR.Libs.thermoFunctions import Sutherland_Law
from NTR.Libs.boundaryLayerFunctions import calcCf
from NTR.PostProcessing.Tecplot.Output.writeTecplotFile import writeTecplot1DFile
from NTR.PostProcessing.Tecplot.Output.writeTecplotFile import writeTecplot2D3DFile
from NTR.Libs.aeroFunctions import Ma
from NTR.Libs.aeroFunctions import p_t_is
from NTR.PostProcessing.OF.Calc.getBladeSurface import getBladeSurface
"""

from scipy import interpolate as itp
from collections import OrderedDict
import numpy as np
import sys
import vtk.util.numpy_support as VN
import matplotlib.pyplot as plt

"""Einlesen des Mittelschnitts"""


def createBoundaryLayerData(case, path_midspan_slice, beta_01, beta_02, midspan_z, Ts, As, incom_bool=False):
    def getProfileCoords(path_slice):

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(path_slice)
        reader.Update()

        polyData = vtk.vtkPolyData()

        appendFilter = vtk.vtkDataSetSurfaceFilter()

        if vtk.VTK_MAJOR_VERSION <= 5:
            appendFilter.SetInput(reader.GetOutput())
        else:
            appendFilter.SetInputData(reader.GetOutput())

        appendFilter.Update()

        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(appendFilter.GetOutput())

        featureEdges = vtk.vtkFeatureEdges()

        if vtk.VTK_MAJOR_VERSION <= 5:
            featureEdges.SetInput(polyData)
        else:
            featureEdges.SetInputData(polyData)

        featureEdges.BoundaryEdgesOn()
        featureEdges.FeatureEdgesOn()
        featureEdges.ManifoldEdgesOff()
        featureEdges.NonManifoldEdgesOff()
        featureEdges.Update()

        data2 = polyData
        data2.BuildLinks()

        data = featureEdges.GetOutput()
        data.BuildLinks()
        bounds = data.GetBounds()

        kDTree = vtk.vtkKdTreePointLocator()
        kDTree.SetDataSet(data)
        kDTree.BuildLocator()
        id = kDTree.FindClosestPoint([bounds[0], bounds[2], bounds[-1]])

        # id=kDTree.FindClosestPoint([bounds[1],bounds[2],bounds[-1]])

        kDTree2 = vtk.vtkKdTreePointLocator()
        kDTree2.SetDataSet(data2)
        kDTree2.BuildLocator()

        idFilter = vtk.vtkIdFilter()

        if vtk.VTK_MAJOR_VERSION <= 5:
            idFilter.SetInput(data)
        else:
            idFilter.SetInputData(data)

        idFilter.SetIdsArrayName("ids")
        idFilter.Update()

        point_ids = []

        pointIds = idFilter.GetOutput().GetPointData().GetArray("ids")
        for i in range(pointIds.GetNumberOfTuples()):
            point_ids.append(int(pointIds.GetTuple(i)[0]))

        outer_IdList = vtk.vtkIdList()
        outer_IdList.InsertNextId(id)

        for i in range(len(point_ids)):

            now_id = int(outer_IdList.GetId(outer_IdList.GetNumberOfIds() - 1))
            now_point = data.GetPoint(now_id)

            for j in range(len(point_ids)):

                now_id_2 = point_ids[j]
                now_point_2 = data.GetPoint(now_id_2)

                if data2.IsEdge(kDTree2.FindClosestPoint([now_point[0], now_point[1], now_point[2]]),
                                kDTree2.FindClosestPoint(
                                    [now_point_2[0], now_point_2[1], now_point_2[2]])) != 0 and outer_IdList.IsId(
                    point_ids[j]) == -1:
                    outer_IdList.InsertNextId(point_ids[j])
                    break

        inner_IdList = vtk.vtkIdList()

        for i in range(len(point_ids)):
            if outer_IdList.IsId(point_ids[i]) == -1:
                inner_IdList.InsertNextId(point_ids[i])

        mapper = vtk.vtkCellDataToPointData()

        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.AddInput(data)
        else:
            mapper.AddInputData(data)

        mapper.Update()

        x = []
        y = []

        xx = []
        yy = []

        for i in range(outer_IdList.GetNumberOfIds()):
            point = data.GetPoint(outer_IdList.GetId(i))
            x.append(point[0])
            y.append(point[1])

        for i in range(inner_IdList.GetNumberOfIds()):
            point = data.GetPoint(inner_IdList.GetId(i))
            xx.append(point[0])
            yy.append(point[1])

        return xx, yy

    def sortPoints(x, y, ind_vk, ind_hk):

        # Punkte der Saugseite bestimmen

        # ersten Punkt nach der Vorderkante bestimmen

        dists = []
        indizes = []

        for i in range(len(x)):
            if i != ind_vk and y[i] > y[ind_vk]:
                dists.append(((x[i] - x[ind_vk]) ** 2 + (y[i] - y[ind_vk]) ** 2) * (0.5))
                indizes.append(i)

        ind_ss_p2 = indizes[dists.index(min(dists))]

        indizes = range(len(x))

        indizes.remove(ind_vk)
        indizes.remove(ind_ss_p2)

        indizes_ss = []

        indizes_ss.append(ind_vk)
        indizes_ss.append(ind_ss_p2)

        ind = ind_ss_p2

        while ind != ind_hk:

            dists = []
            inds = []

            point = (x[ind], y[ind])
            for i in range(len(indizes)):
                point2 = (x[indizes[i]], y[indizes[i]])
                dist = ((point2[0] - point[0]) ** 2 + (point2[1] - point[1]) ** 2) * (0.5)

                if indizes[i] not in indizes_ss:
                    dists.append(dist)
                    inds.append(indizes[i])

            indizes_ss.append(inds[dists.index(min(dists))])
            indizes.remove(inds[dists.index(min(dists))])
            ind = inds[dists.index(min(dists))]

        x_ss = []
        y_ss = []

        for i in range(len(indizes_ss)):
            x_ss.append(x[indizes_ss[i]])
            y_ss.append(y[indizes_ss[i]])

        dists = []
        indizes = []

        for i in range(len(x)):
            if i != ind_vk and y[i] < y[ind_vk]:
                dists.append(((x[i] - x[ind_vk]) ** 2 + (y[i] - y[ind_vk]) ** 2) * (0.5))
                indizes.append(i)

        ind_ps_p2 = indizes[dists.index(min(dists))]

        indizes = range(len(x))

        indizes.remove(ind_ps_p2)

        indizes_ps = []

        indizes_ps.append(ind_vk)
        indizes_ps.append(ind_ps_p2)

        ind = ind_ps_p2

        while ind != ind_hk:

            dists = []
            inds = []

            point = (x[ind], y[ind])
            for i in range(len(indizes)):
                point2 = (x[indizes[i]], y[indizes[i]])
                dist = ((point2[0] - point[0]) ** 2 + (point2[1] - point[1]) ** 2) * (0.5)

                if indizes[i] not in indizes_ps:
                    dists.append(dist)
                    inds.append(indizes[i])

            indizes_ps.append(inds[dists.index(min(dists))])
            indizes.remove(inds[dists.index(min(dists))])
            ind = inds[dists.index(min(dists))]

        x_ps = []
        y_ps = []

        for i in range(len(indizes_ps)):
            x_ps.append(x[indizes_ps[i]])
            y_ps.append(y[indizes_ps[i]])

        return x_ss, y_ss, x_ps, y_ps

    def calcNormalLinesAndPoints(x_ss, y_ss, x_ps, y_ps, dist, nop):

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(path_midspan_slice)
        reader.Update()

        polyData = vtk.vtkPolyData()

        appendFilter = vtk.vtkDataSetSurfaceFilter()
        if vtk.VTK_MAJOR_VERSION <= 5:
            appendFilter.SetInput(reader.GetOutput())
        else:
            appendFilter.SetInputData(reader.GetOutput())
        appendFilter.Update()

        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(appendFilter.GetOutput())

        bounds = polyData.GetBounds()

        kDTree = vtk.vtkKdTreePointLocator()
        kDTree.SetDataSet(polyData)
        kDTree.BuildLocator()

        idFilter = vtk.vtkIdFilter()
        if vtk.VTK_MAJOR_VERSION <= 5:
            idFilter.SetInput(polyData)
        else:
            idFilter.SetInputData(polyData)
        idFilter.SetIdsArrayName("ids")
        idFilter.Update()

        mapper = vtk.vtkCellDataToPointData()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.AddInput(polyData)
        else:
            mapper.AddInputData(polyData)
        mapper.Update()

        #        U_field=mapper.GetOutput().GetPointData().GetArray("U")
        #        V_field=mapper.GetOutput().GetPointData().GetArray("V")
        #        mag_U_field=mapper.GetOutput().GetPointData().GetArray("mag_U")
        #        T_field=mapper.GetOutput().GetPointData().GetArray("T")
        #        rho_field=mapper.GetOutput().GetPointData().GetArray("rho")
        #        p_field=mapper.GetOutput().GetPointData().GetArray("p")

        def calcDeri(x_ss, y_ss):

            tck, u = itp.splprep([x_ss, y_ss], u=None, k=3, s=0.0, per=0)
            num_values = itp.splev(u, tck, der=1)

            return num_values[0], num_values[1]

        def deleteDoublePoints(x, y, z):

            tmp = OrderedDict()
            for point in zip(x, y, z):
                tmp.setdefault(point[:2], point)

            mypoints = tmp.values()

            x_new = []
            y_new = []
            z_new = []

            for i in range(len(mypoints)):
                x_new.append(mypoints[i][0])
                y_new.append(mypoints[i][1])
                z_new.append(mypoints[i][2])

            return x_new, y_new, z_new

        def calcSSValues(x_ss, y_ss):

            xp_num, yp_num = calcDeri(x_ss, y_ss)

            x = []
            y = []
            z = []
            U = []
            V = []
            mag_U = []
            T = []
            rho = []
            bl_thickness = []
            mom_thickness = []
            displ_thickness = []
            H12 = []
            mag_U_inf = []
            tau_w = []  # wandschubspannung
            u_tau = []  # wandschubspannungsgeschwindigkeit
            c_f = []  # reibungsbeiwert
            blade_dist = []

            for i in range(len(x_ss)):

                x_values = []
                y_values = []

                p_tot_values = []

                m = yp_num[i] / xp_num[i]
                m = -1 / m
                a = y_ss[i] - m * x_ss[i]

                x1 = x_ss[i] + 1
                y1 = m * x1 + a

                angle = np.arctan((y1 - y_ss[i]) / (x1 - x_ss[i]))

                if m > 0:

                    if xp_num[i] < 0:
                        x_values.append(x_ss[i] - dist * np.cos(angle))
                    else:
                        x_values.append(x_ss[i] + dist * np.cos(angle))
                    y_values.append(m * x_values[-1] + a)

                if m < 0:

                    if xp_num[i] < 0:
                        x_values.append(x_ss[i] + dist * np.cos(angle))
                    else:
                        x_values.append(x_ss[i] - dist * np.cos(angle))

                    y_values.append(m * x_values[-1] + a)

                line = vtk.vtkLineSource()
                line.SetPoint1(x_ss[i], y_ss[i], bounds[-1])
                line.SetPoint2(x_values[-1], y_values[-1], bounds[-1])
                line.SetResolution(nop)
                line.Update()

                probe = vtk.vtkProbeFilter()
                probe.SetInputConnection(line.GetOutputPort())

                if vtk.VTK_MAJOR_VERSION <= 5:
                    probe.SetSource(mapper.GetOutput())
                else:
                    probe.SetSourceData(mapper.GetOutput())

                probe.Update()

                # Get data from the vtk object (probe) to a numpy array
                mag_U_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("mag_U"))
                U_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("U"))
                V_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("V"))
                T_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("T"))
                rho_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("rho"))
                p_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("p"))
                # num_points = probe.GetOutput().GetNumberOfPoints()

                # daten bereinigen
                n_mag_U_pd = []
                n_U_pd = []
                n_V_pd = []
                n_T_pd = []
                n_rho_pd = []
                n_p_pd = []
                x_nn = []
                y_nn = []
                z_nn = []

                for j in range(len(T_pd)):
                    if T_pd[j] > 0.0:
                        x_p, y_p, z_p = probe.GetOutput().GetPoint(j)

                        x_nn.append(x_p)
                        y_nn.append(y_p)
                        z_nn.append(z_p)

                        n_mag_U_pd.append(mag_U_pd[j])
                        n_U_pd.append(U_pd[j])
                        n_V_pd.append(V_pd[j])
                        n_T_pd.append(T_pd[j])
                        n_rho_pd.append(rho_pd[j])
                        n_p_pd.append(p_pd[j])

                mag_U_pd = n_mag_U_pd
                U_pd = n_U_pd
                V_pd = n_V_pd
                T_pd = n_T_pd
                rho_pd = n_rho_pd
                p_pd = n_p_pd

                """IM nachfolgenden Absatz ist vermutlich ein Fehler"""

                for j in range(len(T_pd)):
                    p_tot_values.append(p_t_is(case.kappa, Ma(mag_U_pd[j], case.kappa, case.R_L, T_pd[j]), p_pd[j]))

                max_p_tot = max(p_tot_values)
                max_mag_U = max(mag_U_pd)

                x_point_new = []
                y_point_new = []
                z_point_new = []
                U_values_new = []
                V_values_new = []
                mag_U_values_new = []
                dist_values = []
                T_values_new = []
                rho_values_new = []

                for j in range(len(mag_U_pd)):

                    # if incom_bool==False:
                    if p_tot_values[j] <= 0.99 * max_p_tot:

                        x_point_new.append(x_nn[j])
                        y_point_new.append(y_nn[j])
                        z_point_new.append(z_nn[j])
                        U_values_new.append(U_pd[j])
                        V_values_new.append(V_pd[j])
                        mag_U_values_new.append(mag_U_pd[j])
                        T_values_new.append(T_pd[j])
                        rho_values_new.append(rho_pd[j])

                    else:

                        break
                    # else:

                if not mag_U_values_new:
                    for j in list(range(len(mag_U_pd))):

                        if mag_U_pd[j] <= 0.99 * max_mag_U:

                            x_point_new.append(x_nn[j])
                            y_point_new.append(y_nn[j])
                            z_point_new.append(z_nn[j])
                            U_values_new.append(U_pd[j])
                            V_values_new.append(V_pd[j])
                            mag_U_values_new.append(mag_U_pd[j])
                            T_values_new.append(T_pd[j])
                            rho_values_new.append(rho_pd[j])

                        else:

                            break

                mag_U_inf.append(mag_U_values_new[-1])

                for j in range(len(x_point_new)):

                    if j == 0:
                        dist_values.append(0)
                    else:
                        dist_values.append(dist_values[-1] + np.sqrt(
                            (x_point_new[j] - x_point_new[j - 1]) ** 2 + (y_point_new[j] - y_point_new[j - 1]) ** 2))

                bl_thickness.append(dist_values[-1])

                x.append(x_point_new)
                y.append(y_point_new)
                z.append(z_point_new)
                U.append(U_values_new)
                V.append(V_values_new)
                blade_dist.append(dist_values)
                mag_U.append(mag_U_values_new)
                T.append(T_values_new)
                rho.append(rho_values_new)
                # berechnung der Impulsverlust und Verdraegungsdicke

                h_values_1 = []
                h_values_2 = []

                wall_bool = False
                wall_bool2 = True

                for j in range(len(mag_U_values_new)):
                    h_values_1.append(1 - (mag_U_values_new[j] / mag_U_inf[-1]))
                    h_values_2.append(
                        (1 - (mag_U_values_new[j] / mag_U_inf[-1])) * (mag_U_values_new[j] / mag_U_inf[-1]))

                    if mag_U_values_new[j] > 0 and wall_bool2 == True:
                        wall_bool = True
                        wall_bool2 = False
                    if wall_bool == True:
                        eta_wand = Sutherland_Law(T_values_new[j], As, Ts)
                        delta_y = np.sqrt(
                            (x_point_new[j + 1] - x_point_new[j]) ** 2 + (y_point_new[j + 1] - y_point_new[j]) ** 2)
                        tau_w.append(eta_wand * (mag_U_values_new[j] / delta_y))
                        u_tau.append(np.sqrt(tau_w[-1] / rho_values_new[j]))
                        c_f.append(tau_w[-1] / (0.5 * rho[-1] * mag_U[-1] ** 2))

                        wall_bool = False

                mom_thickness.append(np.trapz(y=h_values_2, x=dist_values))
                displ_thickness.append(np.trapz(y=h_values_1, x=dist_values))
                H12.append(displ_thickness[-1] / mom_thickness[-1])

            return x, y, z, U, V, mag_U, rho, T, blade_dist, mag_U_inf, mom_thickness, displ_thickness, H12, tau_w, u_tau, c_f

        def calcPSValues(x_ps, y_ps):

            xp_num, yp_num = calcDeri(x_ps, y_ps)

            x = []
            y = []
            z = []
            U = []
            V = []
            mag_U = []
            T = []
            rho = []
            bl_thickness = []
            mom_thickness = []
            displ_thickness = []
            H12 = []
            mag_U_inf = []
            tau_w = []  # wandschubspannung
            u_tau = []  # wandschubspannungsgeschwindigkeit
            c_f = []  # reibungsbeiwert
            blade_dist = []

            for i in range(len(x_ps)):

                x_values = []
                y_values = []

                p_tot_values = []

                m = yp_num[i] / xp_num[i]
                m = -1 / m
                a = y_ps[i] - m * x_ps[i]

                x1 = x_ps[i] + 1
                y1 = m * x1 + a

                angle = np.arctan((y1 - y_ps[i]) / (x1 - x_ps[i]))

                if m > 0:

                    if xp_num[i] >= 0:
                        x_values.append(x_ps[i] - dist * np.cos(angle))
                    else:
                        x_values.append(x_ps[i] + dist * np.cos(angle))
                    y_values.append(m * x_values[-1] + a)

                if m < 0:

                    if xp_num[i] >= 0:
                        x_values.append(x_ps[i] + dist * np.cos(angle))
                    else:
                        x_values.append(x_ps[i] - dist * np.cos(angle))

                    y_values.append(m * x_values[-1] + a)

                line = vtk.vtkLineSource()
                line.SetPoint1(x_ps[i], y_ps[i], bounds[-1])
                line.SetPoint2(x_values[-1], y_values[-1], bounds[-1])
                line.SetResolution(nop)
                line.Update()

                probe = vtk.vtkProbeFilter()
                probe.SetInputConnection(line.GetOutputPort())

                if vtk.VTK_MAJOR_VERSION <= 5:
                    probe.SetSource(mapper.GetOutput())
                else:
                    probe.SetSourceData(mapper.GetOutput())

                probe.Update()

                # Get data from the vtk object (probe) to a numpy array
                mag_U_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("mag_U"))
                U_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("U"))
                V_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("V"))
                T_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("T"))
                rho_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("rho"))
                p_pd = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("p"))
                # num_points = probe.GetOutput().GetNumberOfPoints()

                # daten bereinigen
                n_mag_U_pd = []
                n_U_pd = []
                n_V_pd = []
                n_T_pd = []
                n_rho_pd = []
                n_p_pd = []
                x_nn = []
                y_nn = []
                z_nn = []

                for j in range(len(T_pd)):
                    if T_pd[j] > 0.0:
                        x_p, y_p, z_p = probe.GetOutput().GetPoint(j)

                        x_nn.append(x_p)
                        y_nn.append(y_p)
                        z_nn.append(z_p)

                        n_mag_U_pd.append(mag_U_pd[j])
                        n_U_pd.append(U_pd[j])
                        n_V_pd.append(V_pd[j])
                        n_T_pd.append(T_pd[j])
                        n_rho_pd.append(rho_pd[j])
                        n_p_pd.append(p_pd[j])

                mag_U_pd = n_mag_U_pd
                U_pd = n_U_pd
                V_pd = n_V_pd
                T_pd = n_T_pd
                rho_pd = n_rho_pd
                p_pd = n_p_pd

                for j in range(len(T_pd)):
                    p_tot_values.append(p_t_is(case.kappa, Ma(mag_U_pd[j], case.kappa, case.R_L, T_pd[j]), p_pd[j]))

                max_p_tot = max(p_tot_values)
                max_mag_U = max(mag_U_pd)

                # delta_p=max_p_tot-min_p_tot

                #                for j in range(len(p_tot_values)):
                #                    p_tot_values[j]=p_tot_values[j]-min_p_tot

                x_point_new = []
                y_point_new = []
                z_point_new = []
                U_values_new = []
                V_values_new = []
                mag_U_values_new = []
                dist_values = []
                T_values_new = []
                rho_values_new = []

                for j in range(len(mag_U_pd)):

                    if p_tot_values[j] <= 0.99 * max_p_tot:

                        x_point_new.append(x_nn[j])
                        y_point_new.append(y_nn[j])
                        z_point_new.append(z_nn[j])
                        U_values_new.append(U_pd[j])
                        V_values_new.append(V_pd[j])
                        mag_U_values_new.append(mag_U_pd[j])
                        T_values_new.append(T_pd[j])
                        rho_values_new.append(rho_pd[j])

                    else:

                        break

                if not mag_U_values_new:

                    for j in range(len(mag_U_pd)):
                        if mag_U_pd[j] <= 0.99 * max_mag_U:

                            x_point_new.append(x_nn[j])
                            y_point_new.append(y_nn[j])
                            z_point_new.append(z_nn[j])
                            U_values_new.append(U_pd[j])
                            V_values_new.append(V_pd[j])
                            mag_U_values_new.append(mag_U_pd[j])
                            T_values_new.append(T_pd[j])
                            rho_values_new.append(rho_pd[j])

                        else:

                            break
                        #                if not mag_U_values_new:
                #                    x_point_new.append(x_nn[0])
                #                    y_point_new.append(y_nn[0])
                #                    z_point_new.append(z_nn[0])
                #                    U_values_new.append(U_pd[0])
                #                    V_values_new.append(V_pd[0])
                #                    mag_U_values_new.append(mag_U_pd[0])
                #                    T_values_new.append(T_pd[0])
                #                    rho_values_new.append(rho_pd[0])
                #                    x_point_new.append(x_nn[1])
                #                    y_point_new.append(y_nn[1])
                #                    z_point_new.append(z_nn[1])
                #                    U_values_new.append(U_pd[1])
                #                    V_values_new.append(V_pd[1])
                #                    mag_U_values_new.append(mag_U_pd[1])
                #                    T_values_new.append(T_pd[1])
                #                    rho_values_new.append(rho_pd[1])

                mag_U_inf.append(mag_U_values_new[-1])

                for j in range(len(x_point_new)):

                    if j == 0:
                        dist_values.append(0)
                    else:
                        dist_values.append(dist_values[-1] + np.sqrt(
                            (x_point_new[j] - x_point_new[j - 1]) ** 2 + (y_point_new[j] - y_point_new[j - 1]) ** 2))

                bl_thickness.append(dist_values[-1])

                x.append(x_point_new)
                y.append(y_point_new)
                z.append(z_point_new)
                U.append(U_values_new)
                V.append(V_values_new)
                blade_dist.append(dist_values)
                mag_U.append(mag_U_values_new)
                T.append(T_values_new)
                rho.append(rho_values_new)
                # berechnung der Impulsverlust und Verdraegungsdicke

                h_values_1 = []
                h_values_2 = []

                wall_bool = False
                wall_bool2 = True

                for j in range(len(mag_U_values_new)):
                    h_values_1.append(1 - (mag_U_values_new[j] / mag_U_inf[-1]))
                    h_values_2.append(
                        (1 - (mag_U_values_new[j] / mag_U_inf[-1])) * (mag_U_values_new[j] / mag_U_inf[-1]))

                    if mag_U_values_new[j] > 0 and wall_bool2 == True:
                        wall_bool = True
                        wall_bool2 = False
                    if wall_bool == True:
                        eta_wand = Sutherland_Law(T_values_new[j], As, Ts)
                        delta_y = np.sqrt(
                            (x_point_new[j + 1] - x_point_new[j]) ** 2 + (y_point_new[j + 1] - y_point_new[j]) ** 2)
                        tau_w.append(eta_wand * (mag_U_values_new[j] / delta_y))
                        u_tau.append(np.sqrt(tau_w[-1] / rho_values_new[j]))
                        c_f.append(tau_w[-1] / (0.5 * rho[-1] * mag_U[-1] ** 2))

                        wall_bool = False

                mom_thickness.append(np.trapz(y=h_values_2, x=dist_values))
                displ_thickness.append(np.trapz(y=h_values_1, x=dist_values))
                H12.append(displ_thickness[-1] / mom_thickness[-1])

            return x, y, z, U, V, mag_U, rho, T, blade_dist, mag_U_inf, mom_thickness, displ_thickness, H12, tau_w, u_tau, c_f

        #        def calcSSValues(x_ss, y_ss):
        #
        #            xp_num, yp_num=calcDeri(x_ss, y_ss)
        #
        #            x=[]
        #            y=[]
        #            z=[]
        #            U=[]
        #            V=[]
        #            mag_U=[]
        #            T=[]
        #            rho=[]
        #            bl_thickness=[]
        #            mom_thickness=[]
        #            displ_thickness=[]
        #            H12=[]
        #            mag_U_inf=[]
        #            tau_w=[] #wandschubspannung
        #            u_tau=[] #wandschubspannungsgeschwindigkeit
        #            c_f=[] #reibungsbeiwert
        #            blade_dist=[]
        #
        #            for i in range(len(x_ss)):
        #
        #                x_values=[]
        #                y_values=[]
        #
        #
        #                x_point=[]
        #                y_point=[]
        #                z_point=[]
        #                U_values=[]
        #                V_values=[]
        #                mag_U_values=[]
        #                T_values=[]
        #                rho_values=[]
        #                wall_normal_velo_values=[]Uy
        #                p_tot_values=[]
        #
        #                m=yp_num[i]/xp_num[i]
        #                m=-1/m
        #                a=y_ss[i]-m*x_ss[i]
        #
        #                x1=x_ss[i]+1
        #                y1=m*x1+a
        #
        #                angle=np.arctan((y1-y_ss[i])/(x1-x_ss[i]))
        #
        #                for j in range(int(nop)+1):
        #
        #                    if j>=0:
        #
        #                        if m >0:
        #
        #                            if xp_num[i] <0:
        #                                x_values.append(x_ss[i]-dist*(j/float(nop))*np.cos(angle))
        #                            else:
        #                                x_values.append(x_ss[i]+dist*(j/float(nop))*np.cos(angle))
        #                            y_values.append(m*x_values[-1]+a)
        #
        #                        if m <0:
        #
        #                            if xp_num[i] <0:
        #                                x_values.append(x_ss[i]+dist*(j/float(nop))*np.cos(angle))
        #                            else:
        #                                x_values.append(x_ss[i]-dist*(j/float(nop))*np.cos(angle))
        #
        #                            y_values.append(m*x_values[-1]+a)
        #
        #                for j in range(len(x_values)):
        #
        #                    point_id=kDTree.FindClosestPoint([x_values[j],y_values[j],midspan_z])
        #
        #                    point=polyData.GetPoint(point_id)
        #
        #                    dist_points=np.sqrt((point[0]-x_values[j])**2+(point[1]-y_values[j])**2)
        #
        #                    if dist_points<=0.0002:
        #
        #                        x_point.append(point[0])
        #                        y_point.append(point[1])
        #                        z_point.append(point[2])
        #
        #
        #                x_point,y_point,z_point=deleteDoublePoints(x_point,y_point,z_point)
        #
        #                wall_normal_velo_values=[]
        #
        #                #wall_angle=np.arctan2(y_point[-1]-y_point[0],x_point[-1]-x_point[0])-0.5*np.pi
        #
        #                wall_angle=angle-0.5*np.pi
        #
        #                def rotate(ux,vy, angle):
        #
        #
        #                    ueta = np.cos(angle) * (ux ) - np.sin(angle) * (vy)
        #                    vzeta = np.sin(angle) * (ux) + np.cos(angle) * (vy)
        #                    return ueta , vzeta
        #
        #                for j in range(len(x_point)):
        #
        #                    point_id=kDTree.FindClosestPoint([x_point[j],y_point[j],z_point[j]])
        #                    U_values.append(U_field.GetTuple(point_id)[0])
        #                    V_values.append(V_field.GetTuple(point_id)[0])
        #                    mag_U_values.append(mag_U_field.GetTuple(point_id)[0])
        #                    T_values.append(T_field.GetTuple(point_id)[0])
        #                    rho_values.append(rho_field.GetTuple(point_id)[0])
        #                    wall_normal_velo_values.append(np.fabs(rotate(U_values[-1],V_values[-1], wall_angle)[0]))
        #                    p_tot_values.append(p_t_is(case.kappa,Ma(mag_U_values[-1],case.kappa,case.R_L,T_values[-1]),p_field.GetTuple(point_id)[0]) )
        #                #max_U=max(np.fabs(wall_normal_velo_values))
        #
        #                max_p_tot=max(p_tot_values)
        #                #max_U=max(mag_U_values)
        #
        #                x_point_new=[]
        #                y_point_new=[]
        #                z_point_new=[]
        #                U_values_new=[]
        #                V_values_new=[]
        #                mag_U_values_new=[]
        #                dist_values=[]
        #                T_values_new=[]
        #                rho_values_new=[]
        #                wall_normal_velo_values_new=[]
        #
        #                for j in range(len(x_point)):
        #
        #
        #
        #                    if j<2:
        #                        x_point_new.append(x_point[j])
        #                        y_point_new.append(y_point[j])
        #                        z_point_new.append(z_point[j])
        #                        U_values_new.append(U_values[j])
        #                        V_values_new.append(V_values[j])
        #                        mag_U_values_new.append(mag_U_values[j])
        #                        T_values_new.append(T_values[j])
        #                        rho_values_new.append(rho_values[j])
        #                        wall_normal_velo_values_new.append(wall_normal_velo_values[j])
        #
        #
        #                    if j>=2:
        #                        if p_tot_values[j] <= 0.99*max_p_tot:
        #                        #if mag_U_values[j] <= 0.99*max_U:
        #                            x_point_new.append(x_point[j])
        #                            y_point_new.append(y_point[j])
        #                            z_point_new.append(z_point[j])
        #                            U_values_new.append(U_values[j])
        #                            V_values_new.append(V_values[j])
        #                            mag_U_values_new.append(mag_U_values[j])
        #                            T_values_new.append(T_values[j])
        #                            rho_values_new.append(rho_values[j])
        #                            index_inf=j
        #
        #                        else:
        #
        #                            break
        #
        #                mag_U_inf.append(mag_U_values[index_inf])
        #
        #                #print(len(x_point_new))
        #
        #                for j in range(len(x_point_new)):
        #
        #                    if j==0:
        #                        dist_values.append(0)
        #                    else:
        #                        dist_values.append(dist_values[-1]+np.sqrt((x_point_new[j]-x_point_new[j-1])**2+(y_point_new[j]-y_point_new[j-1])**2))
        #
        #                #print(len(dist_values))
        #
        #                bl_thickness.append(dist_values[-1])
        #
        #                x.append(x_point_new)
        #                y.append(y_point_new)
        #                z.append(z_point_new)
        #                U.append(U_values_new)
        #                V.append(V_values_new)
        #                blade_dist.append(dist_values)
        #                mag_U.append(mag_U_values_new)
        #                T.append(T_values_new)
        #                rho.append(rho_values_new)
        #                #berechnung der Impulsverlust und Verdraegungsdicke
        #
        #                h_values_1=[]
        #                h_values_2=[]
        #
        #                for j in range(len(mag_U_values_new)):
        #                    h_values_1.append(1-(mag_U_values_new[j]/mag_U_inf[-1]))
        #                    h_values_2.append((1-(mag_U_values_new[j]/mag_U_inf[-1]))*(mag_U_values_new[j]/mag_U_inf[-1]))
        #
        #                    if j==0:
        #                        eta_wand=Sutherland_Law(T_values[j],As,Ts)
        #                        delta_y=np.sqrt((x_point[j+1]-x_point[j])**2+(y_point[j+1]-y_point[j])**2)
        #                        tau_w.append(eta_wand*(mag_U_values_new[j]/delta_y))
        #                        u_tau.append(np.sqrt(tau_w[-1]/rho_values_new[j]))
        #                        c_f.append(calcCf(u_tau[-1],mag_U_inf[-1]))
        #
        #                mom_thickness.append(np.trapz(y=h_values_2,x=dist_values))
        #                displ_thickness.append(np.trapz(y=h_values_1,x=dist_values))
        #                H12.append(displ_thickness[-1]/mom_thickness[-1])
        #
        #
        #            return x,y,z,U,V,mag_U,rho,T,blade_dist,mag_U_inf,mom_thickness,displ_thickness,H12,tau_w,u_tau,c_f

        #        def calcPSValues(x_ps, y_ps):
        #
        #            xp_num, yp_num=calcDeri(x_ps, y_ps)
        #
        #            x=[]
        #            y=[]
        #            z=[]
        #            U=[]
        #            V=[]
        #            mag_U=[]
        #            T=[]
        #            rho=[]
        #            bl_thickneps=[]
        #            mom_thickneps=[]
        #            displ_thickneps=[]
        #            H12=[]
        #            mag_U_inf=[]
        #            tau_w=[] #wandschubspannung
        #            u_tau=[] #wandschubspannungsgeschwindigkeit
        #            c_f=[] #reibungsbeiwert
        #            blade_dist=[]
        #
        #
        #            for i in range(len(x_ps)):
        #
        #                x_values=[]
        #                y_values=[]
        #
        #
        #                x_point=[]
        #                y_point=[]
        #                z_point=[]
        #                U_values=[]
        #                V_values=[]
        #                mag_U_values=[]
        #                T_values=[]
        #                rho_values=[]
        #                wall_normal_velo_values=[]
        #                p_tot_values=[]
        #
        #
        #
        #                m=yp_num[i]/xp_num[i]
        #                m=-1/m
        #                a=y_ps[i]-m*x_ps[i]
        #
        #                x1=x_ps[i]+1
        #                y1=m*x1+a
        #
        #                angle=np.arctan((y1-y_ps[i])/(x1-x_ps[i]))
        #
        #                for j in range(int(nop)+1):
        #
        #                    if j>=0:
        #
        #                        if m >0:
        #
        #                            if xp_num[i] >=0:
        #                                x_values.append(x_ps[i]-dist*(j/float(nop))*np.cos(angle))
        #                            else:
        #                                x_values.append(x_ps[i]+dist*(j/float(nop))*np.cos(angle))
        #                            y_values.append(m*x_values[-1]+a)
        #
        #                        if m <0:
        #
        #                            if xp_num[i] >=0:
        #                                x_values.append(x_ps[i]+dist*(j/float(nop))*np.cos(angle))
        #                            else:
        #                                x_values.append(x_ps[i]-dist*(j/float(nop))*np.cos(angle))
        #
        #                            y_values.append(m*x_values[-1]+a)
        #
        #                for j in range(len(x_values)):
        #
        #                    point_id=kDTree.FindClosestPoint([x_values[j],y_values[j],midspan_z])
        #
        #                    point=polyData.GetPoint(point_id)
        #
        #                    dist_points=np.sqrt((point[0]-x_values[j])**2+(point[1]-y_values[j])**2)
        #
        #                    if dist_points<=0.0002:
        #
        #                        x_point.append(point[0])
        #                        y_point.append(point[1])
        #                        z_point.append(point[2])
        #
        #
        #                x_point,y_point,z_point=deleteDoublePoints(x_point,y_point,z_point)
        #
        #                #wall_angle=np.arctan2(y_point[-1]-y_point[0],x_point[-1]-x_point[0])-0.5*np.pi
        #                wall_angle=angle-0.5*np.pi
        #                #wall_angle=np.arctan((y_point[1]-y_point[0])/(x_point[1]-x_point[0]))
        #
        #                def rotate(ux,vy, angle):
        #
        #
        #                    ueta = np.cos(angle) * (ux ) - np.sin(angle) * (vy)
        #                    vzeta = np.sin(angle) * (ux) + np.cos(angle) * (vy)
        #                    return ueta , vzeta
        #
        #                for j in range(len(x_point)):
        #
        #                    point_id=kDTree.FindClosestPoint([x_point[j],y_point[j],z_point[j]])
        #                    U_values.append(U_field.GetTuple(point_id)[0])
        #                    V_values.append(V_field.GetTuple(point_id)[0])
        #                    mag_U_values.append(mag_U_field.GetTuple(point_id)[0])
        #                    T_values.append(T_field.GetTuple(point_id)[0])
        #                    rho_values.append(rho_field.GetTuple(point_id)[0])
        #                    wall_normal_velo_values.append(np.fabs(rotate(U_values[-1],V_values[-1], wall_angle)[0]))
        #                    p_tot_values.append(p_t_is(case.kappa,Ma(mag_U_values[-1],case.kappa,case.R_L,T_values[-1]),p_field.GetTuple(point_id)[0]) )
        #                #max_U=max(np.fabs(wall_normal_velo_values))
        #
        #                #max_U=max(mag_U_values)
        #                max_p_tot=max(p_tot_values)
        #
        #                x_point_new=[]
        #                y_point_new=[]
        #                z_point_new=[]
        #                U_values_new=[]
        #                V_values_new=[]
        #                mag_U_values_new=[]
        #                dist_values=[]
        #                T_values_new=[]
        #                rho_values_new=[]
        #                wall_normal_velo_values_new=[]
        #
        #
        #
        #                for j in range(len(x_point)):
        #
        #
        #
        #                    if j<2:
        #                        x_point_new.append(x_point[j])
        #                        y_point_new.append(y_point[j])
        #                        z_point_new.append(z_point[j])
        #                        U_values_new.append(U_values[j])
        #                        V_values_new.append(V_values[j])
        #                        mag_U_values_new.append(mag_U_values[j])
        #                        T_values_new.append(T_values[j])
        #                        rho_values_new.append(rho_values[j])
        #                        wall_normal_velo_values_new.append(wall_normal_velo_values[j])
        #
        #
        #                    if j>=2:
        #
        #
        #                        #if mag_U_values[j] <= 0.99*max_U:
        #                        if p_tot_values[j] <= 0.99*max_p_tot:
        #                            x_point_new.append(x_point[j])
        #                            y_point_new.append(y_point[j])
        #                            z_point_new.append(z_point[j])
        #                            U_values_new.append(U_values[j])
        #                            V_values_new.append(V_values[j])
        #                            mag_U_values_new.append(mag_U_values[j])
        #                            T_values_new.append(T_values[j])
        #                            rho_values_new.append(rho_values[j])
        #                            index_inf=j
        #
        #                        else:
        #
        #                            break
        #
        #                mag_U_inf.append(mag_U_values[index_inf])
        #
        #                for j in range(len(x_point_new)):
        #                    if j==0:
        #                        dist_values.append(0)
        #                    else:
        #                        dist_values.append(dist_values[-1]+np.sqrt((x_point_new[j]-x_point_new[j-1])**2+(y_point_new[j]-y_point_new[j-1])**2))
        #
        #                bl_thickneps.append(dist_values[-1])
        #
        #                x.append(x_point_new)
        #                y.append(y_point_new)
        #                z.append(z_point_new)
        #                U.append(U_values_new)
        #                V.append(V_values_new)
        #                blade_dist.append(dist_values)
        #                mag_U.append(mag_U_values_new)
        #                T.append(T_values_new)
        #                rho.append(rho_values_new)
        #                #berechnung der Impulsverlust und Verdraegungsdicke
        #
        #                h_values_1=[]
        #                h_values_2=[]
        #
        #                for j in range(len(mag_U_values_new)):
        #                    h_values_1.append(1-(mag_U_values_new[j]/mag_U_inf[-1]))
        #                    h_values_2.append((1-(mag_U_values_new[j]/mag_U_inf[-1]))*(mag_U_values_new[j]/mag_U_inf[-1]))
        #
        #                    if j==0:
        #                        eta_wand=Sutherland_Law(T_values[j],As,Ts)
        #                        delta_y=np.sqrt((x_point[j+1]-x_point[j])**2+(y_point[j+1]-y_point[j])**2)
        #                        tau_w.append(eta_wand*(mag_U_values_new[j]/delta_y))
        #                        u_tau.append(np.sqrt(tau_w[-1]/rho_values_new[j]))
        #                        c_f.append(calcCf(u_tau[-1],mag_U_inf[-1]))
        #
        #                mom_thickneps.append(np.trapz(y=h_values_2,x=dist_values))
        #                displ_thickneps.append(np.trapz(y=h_values_1,x=dist_values))
        #                H12.append(displ_thickneps[-1]/mom_thickneps[-1])
        #
        #
        #            return x,y,z,U,V,mag_U,rho,T,blade_dist,mag_U_inf,mom_thickneps,displ_thickneps,H12,tau_w,u_tau,c_f

        # print(mag_U_bl_ss[100])

        # plt.plot(mag_U_bl_ss[100],y_bl_ss[100],'x-b')

        # plt.plot(x_ss,H12_bl_ss,'-b')
        # plt.plot(x_ss[:575],H12_bl_ss[:575],'-r')

        x_bl_ss, y_bl_ss, z_bl_ss, U_bl_ss, V_bl_ss, mag_U_bl_ss, rho_bl_ss, T_bl_ss, blade_dist_bl_ss, mag_U_inf_bl_ss, \
        mom_thickness_bl_ss, displ_thickness_bl_ss, H12_bl_ss, tau_w_bl_ss, u_tau_bl_ss, c_f_bl_ss = calcSSValues(x_ss,
                                                                                                                  y_ss)
        x_bl_ps, y_bl_ps, z_bl_ps, U_bl_ps, V_bl_ps, mag_U_bl_ps, rho_bl_ps, T_bl_ps, blade_dist_bl_ps, mag_U_inf_bl_ps, \
        mom_thickness_bl_ps, displ_thickness_bl_ps, H12_bl_ps, tau_w_bl_ps, u_tau_bl_ps, c_f_bl_ps = calcPSValues(x_ps,
                                                                                                                  y_ps)

        return x_bl_ss, y_bl_ss, z_bl_ss, U_bl_ss, V_bl_ss, mag_U_bl_ss, rho_bl_ss, T_bl_ss, blade_dist_bl_ss, mag_U_inf_bl_ss, \
               mom_thickness_bl_ss, displ_thickness_bl_ss, H12_bl_ss, tau_w_bl_ss, u_tau_bl_ss, c_f_bl_ss, \
               x_bl_ps, y_bl_ps, z_bl_ps, U_bl_ps, V_bl_ps, mag_U_bl_ps, rho_bl_ps, T_bl_ps, blade_dist_bl_ps, mag_U_inf_bl_ps, \
               mom_thickness_bl_ps, displ_thickness_bl_ps, H12_bl_ps, tau_w_bl_ps, u_tau_bl_ps, c_f_bl_ps

    x_koords, y_koords = getProfileCoords(path_midspan_slice)
    ind_vk, ind_hk = calc_vk_hk(x_koords, y_koords, beta_01, beta_02)
    x_ss, y_ss, x_ps, y_ps = sortPoints(x_koords, y_koords, ind_vk, ind_hk)

    x_bl_ss, y_bl_ss, z_bl_ss, U_bl_ss, V_bl_ss, mag_U_bl_ss, rho_bl_ss, T_bl_ss, blade_dist_bl_ss, mag_U_inf_bl_ss, \
    mom_thickness_bl_ss, displ_thickness_bl_ss, H12_bl_ss, tau_w_bl_ss, u_tau_bl_ss, c_f_bl_ss, \
    x_bl_ps, y_bl_ps, z_bl_ps, U_bl_ps, V_bl_ps, mag_U_bl_ps, rho_bl_ps, T_bl_ps, blade_dist_bl_ps, mag_U_inf_bl_ps, \
    mom_thickness_bl_ps, displ_thickness_bl_ps, H12_bl_ps, tau_w_bl_ps, u_tau_bl_ps, c_f_bl_ps = calcNormalLinesAndPoints(
        x_ss, y_ss, x_ps, y_ps, 0.5 * (x_ss[-1] - x_ss[0]), 8000)

    x_zu_l_ax_ss = []
    x_zu_l_ax_ps = []

    for i in range(len(x_ss)):
        x_zu_l_ax_ss.append((x_ss[i] - x_ss[0]) / (x_ss[-1] - x_ss[0]))

    for i in range(len(x_ps)):
        x_zu_l_ax_ps.append((x_ps[i] - x_ps[0]) / (x_ps[-1] - x_ps[0]))

    def writeOutput():

        output_path = os.path.join(os.path.dirname(os.path.abspath(path_midspan_slice)), 'bl_data_int.dat')
        values = [[x_ss, y_ss, x_zu_l_ax_ss, mag_U_inf_bl_ss, mom_thickness_bl_ss, displ_thickness_bl_ss, H12_bl_ss,
                   tau_w_bl_ss, u_tau_bl_ss, c_f_bl_ss],
                  [x_ps, y_ps, x_zu_l_ax_ps, mag_U_inf_bl_ps, mom_thickness_bl_ps, displ_thickness_bl_ps, H12_bl_ps,
                   tau_w_bl_ps, u_tau_bl_ps, c_f_bl_ps]]
        writeTecplot1DFile(output_path,
                           ['X', 'Y', 'x<sub>Ax</sub> / l<sub>Ax</sub>', 'U_inf', 'Theta', 'delta', 'H12', 'tau_w',
                            'u_tau', 'c_f'], ['Saugseite', 'Druckseite'], values, 'Grenzschichtdaten')

        output_path = os.path.join(os.path.dirname(os.path.abspath(path_midspan_slice)), 'bl_data.dat')

        nops = 0

        for i in range(len(x_bl_ss)):
            nops = nops + len(x_bl_ss[i])

        data = open(output_path, 'w')
        data.write('TITLE = "komplette Grenzschichtdaten fuer z=' + str(midspan_z) + '"\n')
        data.write('Variables="X","Y","U","V","mag_U","rho","T","dist_blade"\n')
        data.write('ZONE T="Saugseite" ,I=' + str(int(nops)) + ', F=POINT\n')

        for i in range(len(x_bl_ss)):
            for j in range(len(x_bl_ss[i])):
                data.write(str(x_bl_ss[i][j]) + "\t" + str(y_bl_ss[i][j]) + "\t" + str(U_bl_ss[i][j]) + "\t" + str(
                    V_bl_ss[i][j]) + "\t" + str(mag_U_bl_ss[i][j]) + "\t" + str(rho_bl_ss[i][j]) + "\t" + str(
                    T_bl_ss[i][j]) + "\t" + str(blade_dist_bl_ss[i][j]) + "\n")

        nops = 0

        for i in range(len(x_bl_ps)):
            nops = nops + len(x_bl_ps[i])

        data.write('ZONE T="Druckseite" ,I=' + str(int(nops)) + ', F=POINT\n')

        for i in range(len(x_bl_ps)):
            for j in range(len(x_bl_ps[i])):
                data.write(str(x_bl_ps[i][j]) + "\t" + str(y_bl_ps[i][j]) + "\t" + str(U_bl_ps[i][j]) + "\t" + str(
                    V_bl_ps[i][j]) + "\t" + str(mag_U_bl_ps[i][j]) + "\t" + str(rho_bl_ps[i][j]) + "\t" + str(
                    T_bl_ps[i][j]) + "\t" + str(blade_dist_bl_ps[i][j]) + "\n")

        data.close()

    writeOutput()

    def createControlPlots():
        # position points
        output_path = os.path.join(os.path.dirname(os.path.abspath(path_midspan_slice)),
                                   'kontrollplot_points_position.pdf')

        plt.figure(figsize=(9, 9))

        plt.plot(x_ss, y_ss, '-k')
        plt.plot(x_ps, y_ps, '-r')

        for i in range(len(x_bl_ss)):
            plt.plot(x_bl_ss[i], y_bl_ss[i], 'x-r', ms=0.01, lw=0.5)

        for i in range(len(x_bl_ps)):
            plt.plot(x_bl_ps[i], y_bl_ps[i], 'x-b', ms=0.01, lw=0.5)

        plt.grid()
        plt.axis('equal')
        plt.savefig(output_path)
        plt.close('all')

        # H12
        # saugseite
        output_path = os.path.join(os.path.dirname(os.path.abspath(path_midspan_slice)), 'kontrollplot_H12_ss.pdf')
        plt.figure(figsize=(9, 9))

        plt.plot(x_zu_l_ax_ss, H12_bl_ss, '-k')

        plt.grid()
        plt.savefig(output_path)
        plt.close('all')

        output_path = os.path.join(os.path.dirname(os.path.abspath(path_midspan_slice)), 'kontrollplot_H12_ps.pdf')
        plt.figure(figsize=(9, 9))

        plt.plot(x_zu_l_ax_ps, H12_bl_ps, '-k')

        plt.grid()
        plt.savefig(output_path)
        plt.close('all')

    createControlPlots()


def createBoundaryLayerData2(path_midspan_slice, path_blade_surface, kappa=1.4, R_L=287.058, As=1.458e-06, Ts=110.4,
                             tolerance=0.01):  # ,beta_01,beta_02,midspan_z,Ts,As,incom_bool=False):

    x, y, xx, yy = getGeom2DVTUSLice(path_midspan_slice)
    x_ss, y_ss, x_ps, y_ps = sortProfilPoints(xx, yy)

    # blade surface einlesen
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path_blade_surface)
    reader.Update()
    blade_surface = reader.GetOutput()

    # mittelschnitt einlesen
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path_midspan_slice)
    reader.Update()
    midspan_slice = reader.GetOutput()
    bounds = midspan_slice.GetBounds()

    # pos 1 werte einlesen
    sys.path.append(os.path.dirname(path_midspan_slice))
    from postSlicesValues import get_postSlicesValues

    class gV():

        b = 0

    gV = gV()
    get_postSlicesValues(gV)

    # bestimme normalenvektor
    points = vtk.vtkPoints()
    pointPolyData = vtk.vtkPolyData()
    norm_vecs = []

    for i in range(len(x_ss)):
        points.InsertNextPoint(x_ss[i], y_ss[i], bounds[-1])

    pointPolyData.SetPoints(points)

    probe = vtk.vtkProbeFilter()

    if vtk.VTK_MAJOR_VERSION <= 5:

        probe.SetInput(pointPolyData)
        probe.SetSource(blade_surface)
    else:
        probe.SetInputData(pointPolyData)
        probe.SetSourceData(blade_surface)

    probe.Update()
    probe_data = probe.GetOutput()

    norms_array = probe_data.GetPointData().GetArray('Normals')

    for i in range(probe_data.GetNumberOfPoints()):
        norm_vec = norms_array.GetTuple(i)

        norm_vec = np.array([float(norm_vec[0]), float(norm_vec[1]), float(norm_vec[2])])

        # face normal wird nachfolgend umgedreht und zeigt jetzt von der profiloberflacehe in das stroemungsgebiet
        norm_vec = -norm_vec
        norm_vecs.append(list(norm_vec))

    values_bl_int = []
    points = probe_data.GetPoints()

    # Gehe Saugseite durch
    # gehe jeden punkt durch und schneide 2d schnitt entlang des normalen vektors

    zone_values_ss = []
    zone_name_strings_ss = []

    x_ax_l_ax = []
    mag_U_inf_bl = []
    mom_thickness_bl = []
    displ_thickness_bl = []
    H12_bl = []
    tau_w_bl = []
    u_tau_bl = []
    c_f_bl = []

    x_ss_new = []
    y_ss_new = []

    for i in range(len(x_ss)):

        cut_plane = vtk.vtkPlane()
        cut_plane.SetOrigin(x_ss[i], y_ss[i], bounds[-1])

        rot_norm_vec = list(rotateTensorVector(np.array(norm_vecs[i]), [0, 0, 1], 90.0))

        cut_plane.SetNormal(rot_norm_vec[0], rot_norm_vec[1], rot_norm_vec[2])

        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(cut_plane)
        cutter.SetInputConnection(reader.GetOutputPort())
        cutter.Update()
        data = cutter.GetOutput()

        mapper = vtk.vtkCellDataToPointData()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.AddInput(data)
        else:
            mapper.AddInputData(data)
        mapper.Update()

        cut_data = mapper.GetOutput()

        pt = []
        p = []
        u = []
        v = []
        w = []
        rho = []
        mag_U = []
        T = []
        wall_dist = []
        dudx = []
        dudy = []
        dudz = []
        dvdx = []
        dvdy = []
        dvdz = []
        dwdx = []
        dwdy = []
        dwdz = []
        x = []
        y = []
        z = []

        # print(cut_data.GetNumberOfPoints())

        for j in xrange(cut_data.GetNumberOfPoints()):
            point = cut_data.GetPoint(j)
            delta_origin_x = point[0] - x_ss[i]
            delta_origin_y = point[1] - y_ss[i]

            get_values_bool = False

            if np.sign(delta_origin_x) == np.sign(norm_vecs[i][0]) and np.sign(delta_origin_y) == np.sign(
                norm_vecs[i][1]):
                get_values_bool = True

            dist_point = np.sqrt((point[0] - x_ss[i]) ** 2 + (point[1] - y_ss[i]) ** 2)

            if get_values_bool == True and dist_point not in wall_dist and dist_point <= tolerance:  # and dist_point<=tolerance:
                p.append(float(cut_data.GetPointData().GetArray("p").GetTuple(j)[0]))
                u.append(float(cut_data.GetPointData().GetArray("U").GetTuple(j)[0]))
                v.append(float(cut_data.GetPointData().GetArray("V").GetTuple(j)[0]))
                w.append(float(cut_data.GetPointData().GetArray("W").GetTuple(j)[0]))
                rho.append(float(cut_data.GetPointData().GetArray("rho").GetTuple(j)[0]))
                mag_U.append(float(cut_data.GetPointData().GetArray("mag_U").GetTuple(j)[0]))
                T.append(float(cut_data.GetPointData().GetArray("T").GetTuple(j)[0]))
                pt.append(p_t_is(kappa, Ma(mag_U[-1], kappa, R_L, T[-1]), p[-1]))
                wall_dist.append(np.sqrt((point[0] - x_ss[i]) ** 2 + (point[1] - y_ss[i]) ** 2))
                x.append(point[0])
                y.append(point[1])
                z.append(point[2])
                dudx.append(float(cut_data.GetPointData().GetArray("dudx").GetTuple(j)[0]))
                dudy.append(float(cut_data.GetPointData().GetArray("dudy").GetTuple(j)[0]))
                dudz.append(float(cut_data.GetPointData().GetArray("dudz").GetTuple(j)[0]))
                dvdx.append(float(cut_data.GetPointData().GetArray("dvdx").GetTuple(j)[0]))
                dvdy.append(float(cut_data.GetPointData().GetArray("dvdy").GetTuple(j)[0]))
                dvdz.append(float(cut_data.GetPointData().GetArray("dvdz").GetTuple(j)[0]))
                dwdx.append(float(cut_data.GetPointData().GetArray("dwdx").GetTuple(j)[0]))
                dwdy.append(float(cut_data.GetPointData().GetArray("dwdy").GetTuple(j)[0]))
                dwdz.append(float(cut_data.GetPointData().GetArray("dwdz").GetTuple(j)[0]))
                # werte nach wall_dist sortieren

        # print(z)

        ind_sorted = sorted(range(len(wall_dist)), key=wall_dist.__getitem__)

        p = list(np.array(p)[ind_sorted])
        u = list(np.array(u)[ind_sorted])
        v = list(np.array(v)[ind_sorted])
        w = list(np.array(w)[ind_sorted])
        rho = list(np.array(rho)[ind_sorted])
        mag_U = list(np.array(mag_U)[ind_sorted])
        T = list(np.array(T)[ind_sorted])
        pt = list(np.array(pt)[ind_sorted])
        wall_dist = list(np.array(wall_dist)[ind_sorted])

        dudx = list(np.array(dudx)[ind_sorted])
        dudy = list(np.array(dudy)[ind_sorted])
        dudz = list(np.array(dudz)[ind_sorted])
        dvdx = list(np.array(dvdx)[ind_sorted])
        dvdy = list(np.array(dvdy)[ind_sorted])
        dvdz = list(np.array(dvdz)[ind_sorted])
        dwdx = list(np.array(dwdx)[ind_sorted])
        dwdy = list(np.array(dwdy)[ind_sorted])
        dwdz = list(np.array(dwdz)[ind_sorted])

        x = list(np.array(x)[ind_sorted])
        y = list(np.array(y)[ind_sorted])
        z = list(np.array(z)[ind_sorted])

        # sehr nahe punkte zueinander entfernen
        new_indizies = []

        for j in range(len(wall_dist)):

            if j == 0 or j == 1:
                new_indizies.append(j)
            else:
                # print('yes')
                if (wall_dist[j] - wall_dist[j - 1]) > 0.000000001:
                    # print('yes')
                    new_indizies.append(j)

        ind_sorted = list(new_indizies)
        p = list(np.array(p)[ind_sorted])
        u = list(np.array(u)[ind_sorted])
        v = list(np.array(v)[ind_sorted])
        w = list(np.array(w)[ind_sorted])
        rho = list(np.array(rho)[ind_sorted])
        mag_U = list(np.array(mag_U)[ind_sorted])
        T = list(np.array(T)[ind_sorted])
        pt = list(np.array(pt)[ind_sorted])
        wall_dist = list(np.array(wall_dist)[ind_sorted])

        dudx = list(np.array(dudx)[ind_sorted])
        dudy = list(np.array(dudy)[ind_sorted])
        dudz = list(np.array(dudz)[ind_sorted])
        dvdx = list(np.array(dvdx)[ind_sorted])
        dvdy = list(np.array(dvdy)[ind_sorted])
        dvdz = list(np.array(dvdz)[ind_sorted])
        dwdx = list(np.array(dwdx)[ind_sorted])
        dwdy = list(np.array(dwdy)[ind_sorted])
        dwdz = list(np.array(dwdz)[ind_sorted])

        x = list(np.array(x)[ind_sorted])
        y = list(np.array(y)[ind_sorted])
        z = list(np.array(z)[ind_sorted])

        # uberspringe axiale position wenn keine werte vorhanden

        if len(p) <= 3:
            continue

        x_ss_new.append(x_ss[i])
        y_ss_new.append(y_ss[i])
        x_ax_l_ax.append((x_ss[i] - min(x_ss)) / (max(x_ss) - min(x_ss)))

        delta_pt = max(pt) - p[0]

        for j in range(len(pt)):
            if pt[j] >= p[0] + 0.99 * delta_pt:
                ind_bl_end = j
                break

        #            else:
        #                break

        # legidlich bl data extrahieren:

        p = p[:ind_bl_end + 1]
        u = u[:ind_bl_end + 1]
        v = v[:ind_bl_end + 1]
        w = w[:ind_bl_end + 1]
        rho = rho[:ind_bl_end + 1]
        mag_U = mag_U[:ind_bl_end + 1]
        T = T[:ind_bl_end + 1]
        pt = pt[:ind_bl_end + 1]
        wall_dist = wall_dist[:ind_bl_end + 1]

        dudx = dudx[:ind_bl_end + 1]
        dudy = dudy[:ind_bl_end + 1]
        dudz = dudz[:ind_bl_end + 1]
        dvdx = dvdx[:ind_bl_end + 1]
        dvdy = dvdy[:ind_bl_end + 1]
        dvdz = dvdz[:ind_bl_end + 1]
        dwdx = dwdx[:ind_bl_end + 1]
        dwdy = dwdy[:ind_bl_end + 1]
        dwdz = dwdz[:ind_bl_end + 1]

        x = x[:ind_bl_end + 1]
        y = y[:ind_bl_end + 1]
        z = z[:ind_bl_end + 1]

        # wandgroessen und grenzschichtgroessen bestimmen

        nu = Sutherland_Law(T[0], As, Ts)

        wall_shear_stress_vec = calcWallShearStress(dudx[0], dudy[0], dudz[0], dvdx[0], dvdy[0], dvdz[0], dwdx[0],
                                                    dwdy[0], dwdz[0], np.array(norm_vecs[i]), rho[0], nu, p[0])

        wall_shear_stress = np.linalg.norm(wall_shear_stress_vec)

        if wall_shear_stress_vec[0] < 0:
            wall_shear_stress = -wall_shear_stress

        friction_velocity = (np.abs(wall_shear_stress) / rho[0]) ** (0.5)

        c_f = calcCf(wall_shear_stress, gV.Mag_U_1, gV.rho_1)

        # berechnung der Impulsverlust und Verdraegungsdicke

        h_values_1 = []
        h_values_2 = []

        for j in range(len(mag_U)):
            h_values_1.append(1 - ((mag_U[j] * rho[j]) / (mag_U[-1] * rho[-1])))
            h_values_2.append((1 - (mag_U[j] / mag_U[-1])) * ((mag_U[j] * rho[j]) / (mag_U[-1] * rho[-1])))

        mom_thickneps = np.trapz(y=h_values_2, x=wall_dist)
        displ_thickneps = np.trapz(y=h_values_1, x=wall_dist)

        if np.isnan(displ_thickneps / mom_thickneps):
            H12 = 0

        elif mom_thickneps == 0:
            H12 = 0

        else:
            H12 = displ_thickneps / mom_thickneps

        mag_U_inf_bl.append(mag_U[-1])
        mom_thickness_bl.append(mom_thickneps)
        displ_thickness_bl.append(displ_thickneps)
        H12_bl.append(H12)
        tau_w_bl.append(wall_shear_stress)
        u_tau_bl.append(friction_velocity)
        c_f_bl.append(c_f)

        zone_values_ss.append(
            [x, y, z, pt, p, u, v, w, rho, mag_U, T, wall_dist, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz])
        zone_name_strings_ss.append('Saugseite, x<sub>Ax</sub> / l<sub>Ax</sub> ' + str(x_ax_l_ax[-1]))

    values_bl_int.append(
        [x_ss_new, y_ss_new, x_ax_l_ax, mag_U_inf_bl, mom_thickness_bl, displ_thickness_bl, H12_bl, tau_w_bl, u_tau_bl,
         c_f_bl])

    # Gehe Druckseite durch
    # gehe jeden punkt durch und schneide 2d schnitt entlang des normalen vektors

    zone_values_ps = []
    zone_name_strings_ps = []

    x_ax_l_ax = []
    mag_U_inf_bl = []
    mom_thickness_bl = []
    displ_thickness_bl = []
    H12_bl = []
    tau_w_bl = []
    u_tau_bl = []
    c_f_bl = []

    x_ps_new = []
    y_ps_new = []

    for i in range(len(x_ps)):

        cut_plane = vtk.vtkPlane()
        cut_plane.SetOrigin(x_ps[i], y_ps[i], bounds[-1])

        rot_norm_vec = list(rotateTensorVector(np.array(norm_vecs[i]), [0, 0, 1], 90.0))

        cut_plane.SetNormal(rot_norm_vec[0], rot_norm_vec[1], rot_norm_vec[2])

        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(cut_plane)
        cutter.SetInputConnection(reader.GetOutputPort())
        cutter.Update()
        data = cutter.GetOutput()

        mapper = vtk.vtkCellDataToPointData()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.AddInput(data)
        else:
            mapper.AddInputData(data)
        mapper.Update()

        cut_data = mapper.GetOutput()

        pt = []
        p = []
        u = []
        v = []
        w = []
        rho = []
        mag_U = []
        T = []
        wall_dist = []
        dudx = []
        dudy = []
        dudz = []
        dvdx = []
        dvdy = []
        dvdz = []
        dwdx = []
        dwdy = []
        dwdz = []
        x = []
        y = []
        z = []

        for j in xrange(cut_data.GetNumberOfPoints()):
            point = cut_data.GetPoint(j)
            delta_origin_x = point[0] - x_ps[i]
            delta_origin_y = point[1] - y_ps[i]

            get_values_bool = False

            if np.sign(delta_origin_x) == np.sign(norm_vecs[i][0]) and np.sign(delta_origin_y) == np.sign(
                norm_vecs[i][1]):
                get_values_bool = True

            dist_point = np.sqrt((point[0] - x_ps[i]) ** 2 + (point[1] - y_ps[i]) ** 2)

            if get_values_bool == True and dist_point not in wall_dist and dist_point <= tolerance:  # and dist_point<=tolerance:
                p.append(float(cut_data.GetPointData().GetArray("p").GetTuple(j)[0]))
                u.append(float(cut_data.GetPointData().GetArray("U").GetTuple(j)[0]))
                v.append(float(cut_data.GetPointData().GetArray("V").GetTuple(j)[0]))
                w.append(float(cut_data.GetPointData().GetArray("W").GetTuple(j)[0]))
                rho.append(float(cut_data.GetPointData().GetArray("rho").GetTuple(j)[0]))
                mag_U.append(float(cut_data.GetPointData().GetArray("mag_U").GetTuple(j)[0]))
                T.append(float(cut_data.GetPointData().GetArray("T").GetTuple(j)[0]))

                if not cut_data.GetPointData().GetArray("pt"):

                    pt.append(p_t_is(kappa, Ma(mag_U[-1], kappa, R_L, T[-1]), p[-1]))

                else:

                    pt.append(float(cut_data.GetPointData().GetArray("pt").GetTuple(j)[0]))

                wall_dist.append(np.sqrt((point[0] - x_ps[i]) ** 2 + (point[1] - y_ps[i]) ** 2))
                x.append(point[0])
                y.append(point[1])
                z.append(point[2])
                dudx.append(float(cut_data.GetPointData().GetArray("dudx").GetTuple(j)[0]))
                dudy.append(float(cut_data.GetPointData().GetArray("dudy").GetTuple(j)[0]))
                dudz.append(float(cut_data.GetPointData().GetArray("dudz").GetTuple(j)[0]))
                dvdx.append(float(cut_data.GetPointData().GetArray("dvdx").GetTuple(j)[0]))
                dvdy.append(float(cut_data.GetPointData().GetArray("dvdy").GetTuple(j)[0]))
                dvdz.append(float(cut_data.GetPointData().GetArray("dvdz").GetTuple(j)[0]))
                dwdx.append(float(cut_data.GetPointData().GetArray("dwdx").GetTuple(j)[0]))
                dwdy.append(float(cut_data.GetPointData().GetArray("dwdy").GetTuple(j)[0]))
                dwdz.append(float(cut_data.GetPointData().GetArray("dwdz").GetTuple(j)[0]))
                # werte nach wall_dist sortieren

        # print(z)

        ind_sorted = sorted(range(len(wall_dist)), key=wall_dist.__getitem__)

        p = list(np.array(p)[ind_sorted])
        u = list(np.array(u)[ind_sorted])
        v = list(np.array(v)[ind_sorted])
        w = list(np.array(w)[ind_sorted])
        rho = list(np.array(rho)[ind_sorted])
        mag_U = list(np.array(mag_U)[ind_sorted])
        T = list(np.array(T)[ind_sorted])
        pt = list(np.array(pt)[ind_sorted])
        wall_dist = list(np.array(wall_dist)[ind_sorted])

        dudx = list(np.array(dudx)[ind_sorted])
        dudy = list(np.array(dudy)[ind_sorted])
        dudz = list(np.array(dudz)[ind_sorted])
        dvdx = list(np.array(dvdx)[ind_sorted])
        dvdy = list(np.array(dvdy)[ind_sorted])
        dvdz = list(np.array(dvdz)[ind_sorted])
        dwdx = list(np.array(dwdx)[ind_sorted])
        dwdy = list(np.array(dwdy)[ind_sorted])
        dwdz = list(np.array(dwdz)[ind_sorted])

        x = list(np.array(x)[ind_sorted])
        y = list(np.array(y)[ind_sorted])
        z = list(np.array(z)[ind_sorted])

        # sehr nahe punkte zueinander entfernen
        new_indizies = []

        for j in range(len(wall_dist)):

            if j == 0 or j == 1:
                new_indizies.append(j)
            else:
                # print('yes')
                if (wall_dist[j] - wall_dist[j - 1]) > 0.000000001:
                    # print('yes')
                    new_indizies.append(j)

        ind_sorted = list(new_indizies)
        p = list(np.array(p)[ind_sorted])
        u = list(np.array(u)[ind_sorted])
        v = list(np.array(v)[ind_sorted])
        w = list(np.array(w)[ind_sorted])
        rho = list(np.array(rho)[ind_sorted])
        mag_U = list(np.array(mag_U)[ind_sorted])
        T = list(np.array(T)[ind_sorted])
        pt = list(np.array(pt)[ind_sorted])
        wall_dist = list(np.array(wall_dist)[ind_sorted])

        dudx = list(np.array(dudx)[ind_sorted])
        dudy = list(np.array(dudy)[ind_sorted])
        dudz = list(np.array(dudz)[ind_sorted])
        dvdx = list(np.array(dvdx)[ind_sorted])
        dvdy = list(np.array(dvdy)[ind_sorted])
        dvdz = list(np.array(dvdz)[ind_sorted])
        dwdx = list(np.array(dwdx)[ind_sorted])
        dwdy = list(np.array(dwdy)[ind_sorted])
        dwdz = list(np.array(dwdz)[ind_sorted])

        x = list(np.array(x)[ind_sorted])
        y = list(np.array(y)[ind_sorted])
        z = list(np.array(z)[ind_sorted])

        # uberspringe axiale position wenn keine werte vorhanden

        if len(p) <= 3:
            continue

        x_ps_new.append(x_ps[i])
        y_ps_new.append(y_ps[i])
        x_ax_l_ax.append((x_ss[i] - min(x_ss)) / (max(x_ss) - min(x_ss)))

        delta_pt = max(pt) - p[0]

        for j in range(len(pt)):
            if pt[j] >= p[0] + 0.99 * delta_pt:
                ind_bl_end = j
                break

        #            else:
        #                break

        # legidlich bl data extrahieren:

        p = p[:ind_bl_end + 1]
        u = u[:ind_bl_end + 1]
        v = v[:ind_bl_end + 1]
        w = w[:ind_bl_end + 1]
        rho = rho[:ind_bl_end + 1]
        mag_U = mag_U[:ind_bl_end + 1]
        T = T[:ind_bl_end + 1]
        pt = pt[:ind_bl_end + 1]
        wall_dist = wall_dist[:ind_bl_end + 1]

        dudx = dudx[:ind_bl_end + 1]
        dudy = dudy[:ind_bl_end + 1]
        dudz = dudz[:ind_bl_end + 1]
        dvdx = dvdx[:ind_bl_end + 1]
        dvdy = dvdy[:ind_bl_end + 1]
        dvdz = dvdz[:ind_bl_end + 1]
        dwdx = dwdx[:ind_bl_end + 1]
        dwdy = dwdy[:ind_bl_end + 1]
        dwdz = dwdz[:ind_bl_end + 1]

        x = x[:ind_bl_end + 1]
        y = y[:ind_bl_end + 1]
        z = z[:ind_bl_end + 1]

        # wandgroessen und grenzschichtgroessen bestimmen

        nu = Sutherland_Law(T[0], As, Ts)

        wall_shear_stress_vec = calcWallShearStress(dudx[0], dudy[0], dudz[0], dvdx[0], dvdy[0], dvdz[0], dwdx[0],
                                                    dwdy[0], dwdz[0], np.array(norm_vecs[i]), rho[0], nu, p[0])

        wall_shear_stress = np.linalg.norm(wall_shear_stress_vec)

        if wall_shear_stress_vec[0] < 0:
            wall_shear_stress = -wall_shear_stress

        friction_velocity = (np.abs(wall_shear_stress) / rho[0]) ** (0.5)

        c_f = calcCf(wall_shear_stress, gV.Mag_U_1, gV.rho_1)

        # berechnung der Impulsverlust und Verdraegungsdicke

        h_values_1 = []
        h_values_2 = []

        for j in range(len(mag_U)):
            h_values_1.append(1 - ((mag_U[j] * rho[j]) / (mag_U[-1] * rho[-1])))
            h_values_2.append((1 - (mag_U[j] / mag_U[-1])) * ((mag_U[j] * rho[j]) / (mag_U[-1] * rho[-1])))

        mom_thickneps = np.trapz(y=h_values_2, x=wall_dist)
        displ_thickneps = np.trapz(y=h_values_1, x=wall_dist)

        if np.isnan(displ_thickneps / mom_thickneps):
            H12 = 0

        elif mom_thickneps == 0:
            H12 = 0

        else:
            H12 = displ_thickneps / mom_thickneps

        mag_U_inf_bl.append(mag_U[-1])
        mom_thickness_bl.append(mom_thickneps)
        displ_thickness_bl.append(displ_thickneps)
        H12_bl.append(H12)
        tau_w_bl.append(wall_shear_stress)
        u_tau_bl.append(friction_velocity)
        c_f_bl.append(c_f)

        zone_values_ps.append(
            [x, y, z, pt, p, u, v, w, rho, mag_U, T, wall_dist, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz])
        zone_name_strings_ps.append('Druckseite, x<sub>Ax</sub> / l<sub>Ax</sub> ' + str(x_ax_l_ax[-1]))

    values_bl_int.append(
        [x_ps_new, y_ps_new, x_ax_l_ax, mag_U_inf_bl, mom_thickness_bl, displ_thickness_bl, H12_bl, tau_w_bl, u_tau_bl,
         c_f_bl])

    def writeOutput():

        output_path = os.path.join(os.path.dirname(os.path.abspath(path_midspan_slice)), 'bl_data_int.dat')
        writeTecplot1DFile(output_path,
                           ['X', 'Y', 'x<sub>Ax</sub> / l<sub>Ax</sub>', 'U_inf', 'Theta', 'delta', 'H12', 'tau_w',
                            'u_tau', 'c_f'], ['Saugseite', 'Druckseite'], values_bl_int, 'Grenzschichtdaten')

        output_path = os.path.join(os.path.dirname(os.path.abspath(path_midspan_slice)), 'bl_data.dat')
        writeTecplot1DFile(output_path,
                           ['X', 'Y', 'Z', 'pt', 'p', 'U', 'V', 'W', 'rho', 'mag_U', 'T', 'wall_dist', 'dudx', 'dudy',
                            'dudz', 'dvdx', 'dvdy', 'dvdz', 'dwdx', 'dwdy', 'dwdz'],
                           zone_name_strings_ss + zone_name_strings_ps, zone_values_ss + zone_values_ps,
                           'Grenzschichtdaten')

    writeOutput()

