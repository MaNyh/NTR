# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 14:53:09 2018

@author: ziesse
"""

import os

from NTR.PreProcessing.Profil.openProfile import Profile
from NTR.PreProcessing.Mesh.prepareMesh import Mesh
from NTR.PreProcessing.Mesh.prepareEIZMesh import EIZMesh
from NTR.PreProcessing.Mesh.prepareEIZCascadeMesh import EIZCascadeMesh
from NTR.PreProcessing.Mesh.duplicateCascadeMesh import Mesh as DuplicatedCascadeMesh
from NTR.PostProcessing.Tecplot.Input.openTecplotFile import openTecplotFile

path_script = os.path.dirname(os.path.abspath(__file__))


class Case():

    def __init__(self):

        # Config Werte
        # allgemeine projektangaben
        self.name = None
        self.project_path = None
        # geom para

        self.beta_01 = None
        self.beta_02 = None
        self.t = None

        # aero para
        self.w2_to_w1 = None
        self.Ma_1 = None
        self.T_1 = None
        self.Ma_2th = None
        self.Re_2th = None

        # fluid para
        self.kappa = None
        self.R_L = None
        self.mu = None

        # mesh para
        self.delta_z_in_vk = None
        self.n = None
        self.exp = None
        self.y_plus_values = None
        self.r_g = None
        self.r_n = None

    def createCase(self, project_path):

        os.mkdir(os.path.join(project_path, 'Profil'))

        def createConfig():
            config_file = open(os.path.join(project_path, 'config.py'), 'w')
            config_file.write('#!/usr/bin/env python2\n')
            config_file.write('# -*- coding: utf-8 -*-\n')
            config_file.write('def get_config(case):\n\n')
            config_file.write('\t#allgemeine projektangaben\n')
            config_file.write('\tcase.name=None\n')
            config_file.write('\tcase.project_path=None\n')
            config_file.write('\tcase.midspan_z=None\n')
            config_file.write('\n')
            config_file.write('\t#geom para\n')
            config_file.write('\tcase.beta_01=None\n')
            config_file.write('\tcase.beta_02=None\n')
            config_file.write('\tcase.t=None\n')
            config_file.write('\n')
            config_file.write('\t#aero para\n')
            config_file.write('\tcase.w2_to_w1=None\n')
            config_file.write('\tcase.Ma_1=None\n')
            config_file.write('\tcase.T_1=None\n')
            config_file.write('\tcase.Ma_2th=None\n')
            config_file.write('\tcase.Re_2th=None\n')
            config_file.write('\n')
            config_file.write('\t#fluid para\n')
            config_file.write('\tcase.kappa=1.4\n')
            config_file.write('\tcase.Rs=287.5\n')
            config_file.write('\tcase.mu=None\n')
            config_file.write('\tcase.rho_inf=None\n')
            config_file.write('\tcase.c_p=None\n')
            config_file.write('\n')
            config_file.write('\t#sutherland para\n')
            config_file.write('\tcase.As=1.458e-06\n')
            config_file.write('\tcase.Ts=110.4\n')
            config_file.write('\n')
            config_file.write('\t#mesh para\n')
            config_file.write('\tcase.delta_z_in_vk=None\n')
            config_file.write('\tcase.n=7.0\n')
            config_file.write('\tcase.exp=1.15\n')
            config_file.write('\tcase.y_plus_values=None\n')
            config_file.write('\tcase.r_n=None\n')
            config_file.write('\tcase.delta_y_vk=None\n')
            config_file.write('\tcase.um_rel_vk=None\n')
            config_file.write('\tcase.um_rel_hk=None\n')
            config_file.write('\tcase.delta_i=None\n')
            config_file.write('\tcase.peri_z=None\n')
            config_file.close()

        def createPara():
            para_file = open(os.path.join(project_path, 'para.py'), 'w')
            para_file.write('#!/usr/bin/env python2\n')
            para_file.write('# -*- coding: utf-8 -*-\n')
            para_file.write('def get_para(case):\n\n')
            para_file.close()

        createPara()
        createConfig()

    def getProfileCoords(self, project_path):
        data_path = os.path.join(project_path, 'Profil', 'profile_data.dat')
        data = openTecplotFile(data_path)
        x_ss = data[0][0]
        y_ss = data[0][1]
        x_ps = data[1][0]
        y_ps = data[1][1]
        return x_ss, y_ss, x_ps, y_ps

    def analyzeProfile(self, project_path):
        profile = Profile()
        profile.openProfile(self, project_path)

        para_file = open(os.path.join(project_path, 'para.py'), 'r')
        para_file_lines = para_file.readlines()
        para_file.close()

        bool_data = False

        for i in range(len(para_file_lines)):
            if para_file_lines[i].startswith('\t#profil para'):
                bool_data = True

        if bool_data == False:
            para_file = open(os.path.join(project_path, 'para.py'), 'a')
            para_file.write('\t#profil para\n')
            para_file.write('\tcase.l=' + str(profile.l) + '\n')
            para_file.write('\tcase.l_ax=' + str(profile.l_ax) + '\n')
            para_file.write('\tcase.l_ss=' + str(profile.l_ss) + '\n')
            para_file.write('\tcase.l_ps=' + str(profile.l_ps) + '\n')
            para_file.write('\tcase.p_LE=[' + str(profile.p_LE[0]) + ',' + str(profile.p_LE[1]) + ']\n')
            para_file.write('\tcase.p_TE=[' + str(profile.p_TE[0]) + ',' + str(profile.p_TE[1]) + ']\n')
            para_file.write('\tcase.nop_ss=' + str(len(profile.x_ss)) + '\n')
            para_file.write('\tcase.nop_ps=' + str(len(profile.x_ps)) + '\n')
            para_file.write('\n')
            para_file.write('\t#aero para\n')
            para_file.write('\tcase.w1=' + str(profile.w1) + '\n')
            para_file.write('\tcase.w2=' + str(profile.w2) + '\n')
            para_file.write('\n')
            para_file.write('\t#sim para\n')
            para_file.write('\tcase.t_profile=' + str(profile.t_profile) + '\n')
            para_file.write('\tcase.t_cascade=' + str(profile.t_cascade) + '\n')
            para_file.write('\tcase.t_plus=' + str(profile.t_plus) + '\n')
            para_file.write('\n')

    def createGeometry(self, project_path):
        geom.createGeometry2(self, project_path)

    def createMesh(self, project_path):
        path_createMesh_script = os.path.join(path_script, '..', 'PreProcessing', 'Mesh', 'createMesh.py')
        cmd = "export PPPath='" + project_path + "'" + "\n" + 'igg121 -batch -print -script ' + path_createMesh_script
        os.system(cmd)

    def createMesh2D(self, project_path):
        path_createMesh_script = os.path.join(path_script, '..', 'PreProcessing', 'Mesh', 'createMesh2D.py')
        cmd = "export PPPath='" + project_path + "'" + "\n" + 'igg111 -batch -print -script ' + path_createMesh_script
        os.system(cmd)

    def createEIZMesh(self, project_path):
        path_createMesh_script = os.path.join(path_script, '..', 'PreProcessing', 'Mesh', 'createEIZMesh.py')
        cmd = "export PPPath='" + project_path + "'" + "\n" + 'igg121 -batch -print -script ' + path_createMesh_script
        # p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        os.system(cmd)

    def createEIZCascadeMesh(self, project_path):
        path_createMesh_script = os.path.join(path_script, '..', 'PreProcessing', 'Mesh', 'createEIZCascadeMesh.py')
        cmd = "export PPPath='" + project_path + "'" + "\n" + 'igg121 -batch -print -script ' + path_createMesh_script
        os.system(cmd)

    def calcMeshQualityIGG(self, project_path):
        path_createMesh_script = os.path.join(path_script, '..', 'PreProcessing', 'Mesh', 'calcMeshQualityIGG.py')

        cmd = "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:" + os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                                                 'Libs')
        os.system(cmd)
        cmd = "export PPPath='" + project_path + "'" + "\n" + 'igg111 -batch -print -script ' + path_createMesh_script
        os.system(cmd)

    def calcMeshQualityOF(self, project_path):
        cmd = "checkMesh -case " + os.path.join(project_path, 'Mesh', 'OF') + ' > ' + os.path.join(project_path, 'Mesh',
                                                                                                   'mesh_quali_of.dat')
        os.system(cmd)

    def prepareMesh(self, project_path):
        mesh = Mesh()
        mesh.prepareMesh(self, project_path)

    def duplicateCascadeMesh(self, project_path):
        mesh = DuplicatedCascadeMesh()
        mesh.duplicateCascadeMesh(self, project_path)


    def prepareEIZMesh(self, project_path):
        mesh = EIZMesh()
        mesh.prepareMesh(self, project_path)

    def prepareEIZCascadeMesh(self, project_path):
        mesh = EIZCascadeMesh()
        mesh.prepareMesh(self, project_path)

    def createSimulation(self, project_path, sim_name):

        if not os.path.exists(os.path.join(project_path, 'Simulationen', sim_name)):
            os.makedirs(os.path.join(project_path, 'Simulationen', sim_name))
            os.makedirs(os.path.join(project_path, 'Simulationen', sim_name, '0'))
            os.makedirs(os.path.join(project_path, 'Simulationen', sim_name, 'constant'))
            os.makedirs(os.path.join(project_path, 'Simulationen', sim_name, 'system'))

            config_file = open(os.path.join(project_path, 'Simulationen', sim_name + '_config.py'), 'w')
            config_file.write('#!/usr/bin/env python2\n')
            config_file.write('# -*- coding: utf-8 -*-\n')
            config_file.write('def get_config(sim):\n\n')
            config_file.write('\t#Allgemeine Angaben\n')
            config_file.write('\tsim.type=None\n')
            config_file.write('\t#RB Inlet\n')
            config_file.write('\tsim.U_in_x=None\n')
            config_file.write('\tsim.U_in_y=None\n')
            config_file.write('\tsim.U_in_z=None\n')
            config_file.write('\tsim.R=None\n')
            config_file.write('\tsim.L=None\n\n')
            config_file.write('\t#RB Outlet\n')
            config_file.write('\tsim.p_out=None\n')
            config_file.write('\t"""Auswertungs Vorgaben"""\n')
            config_file.write('\t"stationär"\n')
            config_file.write('\t#Auswertungs Vorgaben\n')

            config_file.close()

        else:
            raise ValueError('Simulation mit diesem Namen bereits vorhanden')


class globalValues():
    def __init__(self):
        pass
