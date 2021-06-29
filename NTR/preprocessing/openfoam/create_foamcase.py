import os

import NTR
from NTR.preprocessing.openfoam.cascadecase_les_filetemplates import les_templates
from NTR.preprocessing.openfoam.filetemplate_utils import get_template_contents
from NTR.utils.functions import yaml_dict_read
from NTR.preprocessing.openfoam.create_probes import create_probe_dicts



def create_cascadecase_les(settings, mainpath, file_templates, probe_templates):

    directories = file_templates.keys()

    casepath = os.path.join(mainpath, "02_Preprocessing")
    templatepath = os.path.abspath(os.path.dirname(NTR.preprocessing.openfoam.cascadecase_les_filetemplates.__file__))

    create_main_directories(casepath, directories)
    create_files(casepath, settings, file_templates, probe_templates, templatepath)


def create_cascadecase_ras(settings, mainpath, file_templates, probe_templates):

    directories = file_templates.keys()

    casepath = os.path.join(mainpath, "02_Preprocessing")
    templatepath = os.path.abspath(os.path.dirname(NTR.preprocessing.openfoam.cascadecase_les_filetemplates.__file__))

    create_main_directories(casepath, directories)
    create_files(casepath, settings, file_templates, probe_templates, templatepath)


def create_foamcase(setting_file):

    settings = yaml_dict_read(setting_file)
    mainpath = os.path.abspath(os.path.dirname(setting_file))

    if settings["case_settings"]["sim_type"] == "openfoam_les":
        file_templates = les_templates.file_templates
        probe_templates = les_templates.probe_templates
        create_cascadecase_les(settings, mainpath, file_templates, probe_templates)
        create_probe_dicts(settings)

    elif settings["case_settings"]["sim_type"] == "openfoam_ras":
        file_templates = les_templates.file_templates
        probe_templates = les_templates.probe_templates
        create_cascadecase_ras(settings, mainpath, file_templates, probe_templates)


def create_files(casepath, settings, files, probe_templates, templatepath):

    probing_settings = settings["probing"]["probes"]
    probes_dict = {}
    for k, v in probing_settings.items():
        if v == True:
            probes_dict[k] = probe_templates[k]

    templates = get_template_contents(templatepath, files)
    for directory, filenames in files.items():
        for file in filenames:
            template_content = templates[directory][file]

            filesettings = settings["case"]["case_parameters"][file]
            if filesettings:
                for key, value in filesettings.items():
                    template_content = template_content.replace("__" + key + "__", value)

            if file == "controlDict":
                for key, value in probes_dict.items():
                    template_content = template_content.replace("//__globalsetting__" + key + "__//", value)
                template_content = template_content.replace("__DELTAT__", str(settings["case_settings"]["timestep"]))
            if file == "createPatchDict":
                template_content = template_content.replace("__ZSPAN__", str(settings["mesh"]["extrudeLength"]))
            if file == "U" and settings["case_settings"]["sim_type"]=="openfoam_les":
                template_content = template_content.replace("__PITCHPER__", str(settings["geometry"]["pitch"]))
                template_content = template_content.replace("__SPANPER__", str(settings["mesh"]["extrudeLength"]))


            with open(os.path.join(casepath, directory, file), "w", newline='\n') as fobj:
                fobj.writelines(template_content)


def create_main_directories(casepath, directories):
    if not os.path.isdir(casepath):
        os.mkdir(casepath)
    for d in directories:
        if not os.path.isdir(os.path.join(casepath, d)):
            os.mkdir(os.path.join(casepath, d))
