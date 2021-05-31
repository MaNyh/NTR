#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 15:14:47 2018

@author: ziesse
"""
import subprocess
import os
import sys
import shutil
import gc


from NTR.PostProcessing.Convert.convertVTKtoTecplotFile import convertVTK3DtoTecplot
from NTR.Libs.simFunctions import sort
from NTR.PostProcessing.OF.Calc.createSlice import createPlanarSlice
from NTR.PostProcessing.OF.Calc.createProfileData import createProfileData
from NTR.PostProcessing.OF.Calc.createBoundaryLayerData import createBoundaryLayerData2
from NTR.PlotCreation.plotVTKSlice import plotVTKSlice2
from NTR.PostProcessing.Convert.convertVTKtoTecplotFile import convertVTKtoVTU
from NTR.Libs.simFunctions import sortOFTimeSteps
from NTR.PostProcessing.OF.Calc.createSpanwiseAverage import createSpanwiseAverage
from NTR.PostProcessing.OF.Calc.getBladeSurface import getBladeSurface
from NTR.PostProcessing.OF.Calc.createNLData import createNLData
from NTR.Libs.geomFunctions import extractCellsByBox



def createStatData(case, sim_path, reconstruct_bool=True, incom_bool=False, foamToVTK_bool=True, bl_bool=True,
                   nl_bool=True, shrink_domain=False):
    incom_bool = incom_bool

    def clean():

        if os.path.exists(os.path.join(sim_path, 'logs', 'run_auswertung.log')):
            os.remove(os.path.join(sim_path, 'logs', 'run_auswertung.log'))

        if os.path.exists(os.path.join(sim_path, 'VTK')):
            shutil.rmtree(os.path.join(sim_path, 'VTK'))

    def reconstruct():

        out = subprocess.check_output(
            "reconstructPar -latestTime -fields '(UMean TMean rhoMean pMean UPrime2Mean nutMean tauMean tauGradUMean gradUMean kMean omegaMean)' -case " + sim_path + " > " + sim_path + '/logs/run_auswertung.log',
            stderr=subprocess.STDOUT, shell=True)

    def calcUMeanGrad():

        out = subprocess.check_output(
            "postProcess -latestTime -func 'grad(UMean)' -case " + sim_path + " > " + sim_path + '/logs/run_auswertung.log',
            stderr=subprocess.STDOUT, shell=True)

    def convertFoamtoVTK():

        out = subprocess.check_output(
            "foamToVTK -latestTime -fields '(UMean TMean rhoMean pMean UPrime2Mean grad(UMean) gradUMean nutMean tauMean tauGradUMean kMean omegaMean)' -case " + sim_path + " > " + sim_path + '/logs/run_auswertung.log',
            stderr=subprocess.STDOUT, shell=True)

    def convertVTKtoTecplot():

        listdir = os.listdir(os.path.join(sim_path, 'VTK'))

        vtk_file_name = ''

        for i in range(len(listdir)):
            if listdir[i].endswith('.vtk'):
                vtk_file_name = listdir[i]

        time_steps = float(vtk_file_name.split('.vtk')[0].split('_')[-1])

        # print(time_steps)

        sim_name = ''

        strings = vtk_file_name.split('.vtk')[0].split('_')

        for i in range(len(strings) - 1):

            if i < len(strings) - 2:
                sim_name = sim_name + strings[i] + '_'
            else:
                sim_name = sim_name + strings[i]

        # print(sim_name)

        t = [i for i in os.listdir(os.path.join(sim_path)) if '0.' in i or 'e-' in i]

        t = sort(t)[-1]

        if not os.path.isdir(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D')):
            os.makedirs(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D'))

        if incom_bool == True:

            convertVTK3DtoTecplot(os.path.join(sim_path, 'VTK', vtk_file_name),
                                  os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                                  case.Ts, case.As, case.kappa, case.R_L, incom_bool=True, rho_incom=case.rho_incom,
                                  T_incom=case.T_incom)

        else:

            convertVTK3DtoTecplot(os.path.join(sim_path, 'VTK', vtk_file_name),
                                  os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                                  case.Ts, case.As, case.kappa, case.R_L, incom_bool=False)

        if shrink_domain == True:
            extractCellsByBox(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                              os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                              case.x_pos_1 - 0.002, 1, -1, 1, -1, 1)

        shutil.rmtree(os.path.join(sim_path, 'VTK'))

    def createMidspanSlice():

        sim_name = sim_path.split('/')[-1]

        t = [i for i in os.listdir(os.path.join(sim_path)) if '0.' in i or 'e-' in i]

        t = sort(t)[-1]

        if not os.path.isdir(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan')):
            os.makedirs(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan'))

        if case.midspan_z == 0:
            cut_z = 0.00001
        else:
            cut_z = case.midspan_z * 0.999

        createPlanarSlice(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                          os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan',
                                       'midspan.vtu'),
                          [0, 0, cut_z], [0, 0, 1])

        createProfileData(case, os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan',
                                             'midspan.vtu'), case.x_pos_1, case.x_pos_2)

        sys.path.append(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan'))

        from postSlicesValues import get_postSlicesValues

        class gV():
            abc = 1

        gV = gV
        get_postSlicesValues(gV)

        getBladeSurface(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                        os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', 'blade_surface.vtu'),
                        gV.Mag_U_1, gV.rho_1, kappa=case.kappa, R_L=case.R_L, tolerance=0.001, Ts=case.Ts, As=case.As)

        if bl_bool == True:
            createBoundaryLayerData2(
                os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan', 'midspan.vtu'),
                os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', 'blade_surface.vtu'),
                kappa=case.kappa, R_L=case.R_L, As=case.As, Ts=case.Ts)

        if nl_bool == True:

            if not os.path.isdir(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan', 'NL')):
                os.makedirs(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan', 'NL'))

            createNLData(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan', 'midspan.vtu'),
                         os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan', 'NL'),
                         case.x_pos_1, case.x_pos_2, case.p_k, case.kappa, case.R_L, case.l, case.As, case.c_p, case.Ts,
                         animation=False, row_number_y_rel=1, row_number_zeta=2)

    #        	plotVTKSlice2(os.path.join('.','Auswertung','Daten',sim_name,'stat',t,'B2B','midspan','midspan.vtu'),
    #        	'mag_U',os.path.join('.','Auswertung','Daten',sim_name,'stat',t,'B2B','midspan','midspan_mag_U.png'),cmap='jet',label_plot='test',label_cb='mag_U',
    #        	num_colorlevels=300,dpi=100,figsize=[10,10])

    def createSpanwiseAverageSlice():

        sim_name = sim_path.split('/')[-1]

        t = [i for i in os.listdir(os.path.join(sim_path)) if '0.' in i or 'e-' in i]

        t = sort(t)[-1]

        if not os.path.isdir(
            os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan_spanwise_average')):
            os.makedirs(
                os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan_spanwise_average'))

        if case.midspan_z == 0:
            cut_z = 0.00001
        else:
            cut_z = case.midspan_z * 0.999

        Origin = [0, 0, cut_z]
        Normal = [0, 0, 1]

        if case.peri_z == True:
            createSpanwiseAverage(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                                  os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B',
                                               'midspan_spanwise_average'),
                                  Origin, Normal, peri_z=case.peri_z)
        else:
            createSpanwiseAverage(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                                  os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B',
                                               'midspan_spanwise_average'),
                                  Origin, Normal, peri_z=case.peri_z, nops=case.nops_spanwise_average)

        createProfileData(case, os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B',
                                             'midspan_spanwise_average', 'midspan_spanwise_average.vtu'), case.x_pos_1,
                          case.x_pos_2)

        if bl_bool == True:
            createBoundaryLayerData2(
                os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan_spanwise_average',
                             'midspan_spanwise_average.vtu'),
                os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', 'blade_surface.vtu'),
                kappa=case.kappa, R_L=case.R_L, As=case.As, Ts=case.Ts)
        if nl_bool == True:

            if not os.path.isdir(
                os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan_spanwise_average', 'NL')):
                os.makedirs(
                    os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan_spanwise_average',
                                 'NL'))

            createNLData(
                os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan_spanwise_average',
                             'midspan_spanwise_average.vtu'),
                os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'B2B', 'midspan_spanwise_average', 'NL'),
                case.x_pos_1, case.x_pos_2, case.p_k, case.kappa, case.R_L, case.l, case.As, case.c_p, case.Ts,
                animation=False, row_number_y_rel=1, row_number_zeta=2)

    def createPos1Slice():

        sim_name = sim_path.split('/')[-1]

        t = [i for i in os.listdir(os.path.join(sim_path)) if '0.' in i or 'e-' in i]

        t = sort(t)[-1]

        if not os.path.isdir(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_1')):
            os.makedirs(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_1'))

        createPlanarSlice(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                          os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_1',
                                       'pos_1_slice.vtu'),
                          [case.x_pos_1, 0, 0], [1, 0, 0])

        plotVTKSlice2(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_1', 'pos_1_slice.vtu'),
                      'mag_U', os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_1',
                                            'pos_1_slice_mag_U.png'), cmap='jet', label_plot='Auswerteebene 1',
                      label_cb='mag_U',
                      num_colorlevels=300, dpi=100, figsize=[10, 10])

    def createPos2Slice():

        sim_name = sim_path.split('/')[-1]

        t = [i for i in os.listdir(os.path.join(sim_path)) if '0.' in i or 'e-' in i]

        t = sort(t)[-1]

        if not os.path.isdir(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_2')):
            os.makedirs(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_2'))

        createPlanarSlice(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, '3D', '3d.vtu'),
                          os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_2',
                                       'pos_2_slice.vtu'),
                          [case.x_pos_2, 0, 0], [1, 0, 0])

        plotVTKSlice2(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_2', 'pos_2_slice.vtu'),
                      'mag_U', os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'POS_2',
                                            'pos_2_slice_mag_U.png'), cmap='jet', label_plot='Auswerteebene 2',
                      label_cb='mag_U',
                      num_colorlevels=300, dpi=100, figsize=[10, 10])

    def getInletSlice():
        sim_name = sim_path.split('/')[-1]

        t = [i for i in os.listdir(os.path.join(sim_path)) if '0.' in i or 'e-' in i]

        sorted_time_steps_float, sorted_time_steps_str = sortOFTimeSteps(t)

        t = sorted_time_steps_str[-1]

        if not os.path.isdir(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'INLET')):
            os.makedirs(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'INLET'))

        boundary_patch_names = subprocess.check_output(
            "foamDictionary -entry boundaryField -keywords " + sim_path + '/' + t + "/UMean -case " + sim_path,
            stderr=subprocess.STDOUT, shell=True)

        boundary_patch_names = boundary_patch_names.rstrip().split()

        # print(boundary_patch_names)

        to_exclude_string = "'("

        for i in range(len(boundary_patch_names)):
            if boundary_patch_names[i] != 'INLET' and boundary_patch_names[i] != 'Inlet' and boundary_patch_names[
                i] != 'inlet':
                to_exclude_string = to_exclude_string + boundary_patch_names[i]

                if i < len(boundary_patch_names) - 1:
                    to_exclude_string = to_exclude_string + ' '

        to_exclude_string = to_exclude_string + ")'"

        # print(sim_path)

        # print("foamToVTK -latestTime -noInternal -excludePatches "+to_exclude_string+" -fields '(UMean TMean rhoMean pMean UPrime2Mean nutMean tauMean tauGradUMean gradUMean kMean omegaMean)' -case "+sim_path)

        out = subprocess.check_output(
            "foamToVTK -latestTime -nearCellValue -noPointValues -noInternal -excludePatches " + to_exclude_string + " -fields '(UMean TMean rhoMean pMean UPrime2Mean nutMean tauMean tauGradUMean gradUMean kMean omegaMean)' -case " + sim_path,
            stderr=subprocess.STDOUT, shell=True)

        vtk_dir_list = os.listdir(os.path.join(sim_path, 'VTK', 'INLET'))

        for i in range(len(vtk_dir_list)):
            if vtk_dir_list[i].startswith('INLET'):
                inlet_vtk_slice_name = vtk_dir_list[i]

        # print(os.path.join(sim_path,'VTK','INLET',inlet_vtk_slice_name))

        convertVTKtoVTU(os.path.join(sim_path, 'VTK', 'INLET', inlet_vtk_slice_name),
                        os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'INLET', 'inlet_slice.vtu'))
        shutil.rmtree(os.path.join(sim_path, 'VTK'))

    def getOutletSlice():
        sim_name = sim_path.split('/')[-1]

        t = [i for i in os.listdir(os.path.join(sim_path)) if '0.' in i or 'e-' in i]

        sorted_time_steps_float, sorted_time_steps_str = sortOFTimeSteps(t)

        t = sorted_time_steps_str[-1]

        if not os.path.isdir(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'OUTLET')):
            os.makedirs(os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'OUTLET'))

        boundary_patch_names = subprocess.check_output(
            "foamDictionary -entry boundaryField -keywords " + sim_path + '/' + t + "/UMean -case " + sim_path,
            stderr=subprocess.STDOUT, shell=True)

        boundary_patch_names = boundary_patch_names.rstrip().split()

        # print(boundary_patch_names)

        to_exclude_string = "'("

        for i in range(len(boundary_patch_names)):
            if boundary_patch_names[i] != 'OUTLET' and boundary_patch_names[i] != 'Outlet' and boundary_patch_names[
                i] != 'outlet':
                to_exclude_string = to_exclude_string + boundary_patch_names[i]

                if i < len(boundary_patch_names) - 1:
                    to_exclude_string = to_exclude_string + ' '

        to_exclude_string = to_exclude_string + ")'"

        # print(sim_path)

        # print("foamToVTK -latestTime -noInternal -excludePatches "+to_exclude_string+" -fields '(UMean TMean rhoMean pMean UPrime2Mean nutMean tauMean tauGradUMean gradUMean kMean omegaMean)' -case "+sim_path)

        out = subprocess.check_output(
            "foamToVTK -latestTime -nearCellValue -noPointValues -noInternal -excludePatches " + to_exclude_string + " -fields '(UMean TMean rhoMean pMean UPrime2Mean nutMean tauMean tauGradUMean gradUMean kMean omegaMean)' -case " + sim_path,
            stderr=subprocess.STDOUT, shell=True)

        vtk_dir_list = os.listdir(os.path.join(sim_path, 'VTK', 'OUTLET'))

        for i in range(len(vtk_dir_list)):
            if vtk_dir_list[i].startswith('OUTLET'):
                inlet_vtk_slice_name = vtk_dir_list[i]

        # print(os.path.join(sim_path,'VTK','OUTLET',inlet_vtk_slice_name))

        convertVTKtoVTU(os.path.join(sim_path, 'VTK', 'OUTLET', inlet_vtk_slice_name),
                        os.path.join('.', 'Auswertung', 'Daten', sim_name, 'stat', t, 'yz', 'OUTLET',
                                     'outlet_slice.vtu'))
        shutil.rmtree(os.path.join(sim_path, 'VTK'))

        # createGlobalData
        # createSlices
        # pos 1 pos 2 slices
        # calc int valus
        # xy_Slices
        ##calc profil data
        ##calc BP data
        ##calc BL data
        ##calc x_lax_data z.B totaldruckverlust
        ##create control plots

        # yz_slices
        # profil_orth_slices

    if foamToVTK_bool == True:
        clean()
        gc.collect()

    if reconstruct_bool == True:
        reconstruct()
        # calcUMeanGrad()
        gc.collect()

    if foamToVTK_bool == True:
        convertFoamtoVTK()
        gc.collect()

    t = convertVTKtoTecplot()
    gc.collect()
    createMidspanSlice()
    gc.collect()
    createSpanwiseAverageSlice()
    gc.collect()
    createPos1Slice()
    gc.collect()
    createPos2Slice()
    gc.collect()

    if reconstruct_bool == True and foamToVTK_bool == True:
        getInletSlice()
        gc.collect()
        getOutletSlice()
        gc.collect()

