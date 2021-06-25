import os

from NTR.preprocessing.openfoam.cascadecase_les_filetemplates.templates import get_template_contents, file_templates, \
    probe_templates


def create_cascadecase_les(settings, mainpath):

    directories = file_templates.keys()

    casepath = os.path.join(mainpath, "02_Preprocessing")

    create_main_directories(casepath, directories)
    create_files(casepath,settings)

def create_case(settings, mainpath):
    if settings[]:
        pass

def create_files(casepath, settings):
    files = file_templates

    probing_settings = settings["probing"]["probes"]
    probes_dict = {}
    for k, v in probing_settings.items():
        if v == True:
            probes_dict[k] = probe_templates[k]

    templates = get_template_contents()
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


            with open(os.path.join(casepath, directory, file), "w", newline='\n') as fobj:
                fobj.writelines(template_content)


def create_main_directories(casepath, directories):
    if not os.path.isdir(casepath):
        os.mkdir(casepath)
    for d in directories:
        if not os.path.isdir(os.path.join(casepath, d)):
            os.mkdir(os.path.join(casepath, d))
