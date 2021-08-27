import os

import NTR
from NTR.preprocessing.case_creation import case_template
from NTR.preprocessing.openfoam.cascadecase_les_filetemplates import les_templates
from NTR.preprocessing.openfoam.create_probes import create_of_les_probe_dicts


class Openfoam_cascade_les(case_template):

    def __init__(self, settings_dict,casepath):
        case_templates = les_templates
        templatepath = os.path.abspath(os.path.dirname(NTR.preprocessing.openfoam.cascadecase_les_filetemplates.__file__))

        super().__init__(case_templates, settings_dict, casepath, templatepath, self.Openfoam_cascade_les_filehandling, create_of_les_probe_dicts)


    def Openfoam_cascade_les_filehandling(self,file, template_content, settings, probes_dict):
        if file == "controlDict":
            for key, value in probes_dict.items():
                template_content = template_content.replace("//__globalsetting__" + key + "__//", value)
            template_content = template_content.replace("__DELTAT__",
                                                        str(settings["openfoam_cascade_les_settings"]["timestep"]))

        if file == "createPatchDict":
            template_content = template_content.replace("__ZSPAN__", str(settings["mesh"]["extrudeLength"]))

        if file == "U":
            template_content = template_content.replace("__PITCHPER__", str(settings["geometry"]["pitch"]))
            template_content = template_content.replace("__SPANPER__", str(settings["mesh"]["extrudeLength"]))
        return template_content


class Openfoam_cascade_ras(case_template):
    def __init__(self, settings_dict,casepath):
        case_templates = les_templates
        templatepath = os.path.abspath(os.path.dirname(NTR.preprocessing.openfoam.cascadecase_ras_filetemplates.__file__))
        super().__init__(case_templates, settings_dict, casepath, templatepath, self.Openfoam_cascade_ras_filehandling, create_of_ras_probe_dicts)
        self.directories = []

    def Openfoam_cascade_ras_filehandling(self,file, template_content, settings, probes_dict):
        return template_content
