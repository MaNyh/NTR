import os

from NTR.utils.functions import yaml_dict_read
from NTR.preprocessing.openfoam.cascadecase_filetemplates.templates import get_template_contents, file_templates


def create_cascadecase_les(yaml_path):

    path = os.path.abspath(os.path.dirname(yaml_path))
    fpath = os.path.join(path,yaml_path)
    settings = yaml_dict_read(fpath)
    directories = file_templates.keys()

    casepath = os.path.join(os.path.dirname(yaml_path), settings["case_name"])


    create_main_directories(casepath, directories)
    create_files(casepath,settings)



def create_files(casepath,settings):
    files = file_templates

    templates = get_template_contents()
    for directory, filenames in files.items():
        for file in filenames:
            template_content = templates[directory][file]
            filesettings = settings["case_parameters"][file]
            if filesettings:
                for key, value in filesettings.items():
                    template_content = template_content.replace("__"+key+"__", value)
            with open(os.path.join(casepath, directory, file), "w") as fobj:
                fobj.writelines(template_content)


def create_main_directories(casepath, directories):
    os.mkdir(casepath)
    for d in directories:
        os.mkdir(os.path.join(casepath, d))


