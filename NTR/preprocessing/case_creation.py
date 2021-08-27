import os

from NTR.preprocessing.openfoam.create_foamcase import Openfoam_cascade_les, Openfoam_cascade_ras
from NTR.utils.filehandling import yaml_dict_read, get_template_contents, read_pickle

from NTR.preprocessing.openfoam.create_probes import create_of_les_probe_dicts



class case_template:
    def __init__(self,case_templates,settings_dict, casepath,templatepath,casespecific_filehandling,create_probe_dicts):
        self.casedirectories = {"ressources":"00_Ressources",
                               "meshing":"01_Meshing",
                               "simcase":"02_Simcase",
                               "solution":"03_Solution",
                                "data":"04_Data"}

        self.settings_dict = settings_dict
        self.casename = settings_dict["case_settings"]["name"]
        self.casepath = casepath
        self.templatepath = templatepath
        self.directories = list(case_templates.file_templates.keys())
        self.filetemplates = case_templates.file_templates
        self.probe_templates = case_templates.probe_templates
        self.sim_directory=os.path.join(self.casepath,self.casedirectories["simcase"])
        self.casespecific_filehandling = casespecific_filehandling
        self.create_probe_dicts = create_probe_dicts


    def create_dirstructure(self):
        for d in self.casedirectories.values():
            if not os.path.isdir(os.path.join(self.casepath,d)):
                os.mkdir(os.path.join(self.casepath,d))
        if not os.path.isdir(self.sim_directory):
            os.mkdir(self.sim_directory)
        for d in self.directories:
            if not os.path.isdir(os.path.join(self.sim_directory, d)):
                os.mkdir(os.path.join(self.sim_directory, d))
        return 0

    def create_files(self):
        probing_settings = self.settings_dict["probing"]["probes"]
        probes_dict = {}
        for k, v in probing_settings.items():
            if v == True:
                probes_dict[k] = self.probe_templates[k]

        templates = get_template_contents(self.templatepath, self.filetemplates)
        for directory, filenames in self.filetemplates.items():
            for file in filenames:
                template_content = templates[directory][file]

                filesettings = self.settings_dict["case"]["case_parameters"][file]
                if filesettings:
                    for key, value in filesettings.items():
                        template_content = template_content.replace("__" + key + "__", value)

                template_content = self.casespecific_filehandling(file,template_content,self.settings_dict,probes_dict)

                with open(os.path.join(self.sim_directory, directory, file), "w", newline='\n') as fobj:
                    fobj.writelines(template_content)

        return 0

    def create_monitors(self, method, args):
        method(*args)
        return 0


def create_simulationcase(path_to_yaml_dict):
    settings = yaml_dict_read(path_to_yaml_dict)
    mainpath = os.path.abspath(os.path.dirname(path_to_yaml_dict))

    case_type = settings["case_settings"]["case_type"]
    geo_ressources = read_pickle(os.path.join(mainpath, "00_Ressources", "01_Geometry", "geometry.pkl"))
    if case_type == "Openfoam_cascade_les":
        casepath = os.path.abspath(os.path.dirname(path_to_yaml_dict))
        Case = Openfoam_cascade_les(settings, casepath)
        Case.create_dirstructure()
        Case.create_files()
        Case.create_monitors(create_of_les_probe_dicts, [settings, geo_ressources])

    if case_type == "Openfoam_cascade_ras":
        Case = Openfoam_cascade_ras(settings, casepath)

        Case.create_dirstructure()
        Case.create_files()
        Case.create_monitors(create_of_les_probe_dicts)
