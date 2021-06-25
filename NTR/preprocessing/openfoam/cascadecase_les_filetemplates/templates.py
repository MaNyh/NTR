import os


def readtxtfile(path_to_file):
    basepath = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(basepath,path_to_file), "r") as fobj:
        content = fobj.readlines()
    return "".join(content)


file_templates = {"0": ["alphat",
                        "nut",
                        "p",
                        "T",
                        "U"],

                  "constant": ["thermophysicalProperties",
                               "turbulenceProperties"],

                  "system": ["controlDict",
                             "createPatchDict",
                             "decomposeParDict",
                             "fvSchemes",
                             "fvSolution",
                             "mapFieldsDict",
                             "Probes_FieldAve_Dict",
                             "Probe_Slice_Dict",
                             "topoSetDict"],

                  "utils": ["monitor.py"],

                  ".": ["submit_job_hlrn.sh"],

                  }

probe_templates = {}

probe_templates["massflow_probing"] = """
    #includeFunc  flowRatePatch(name=INLET)
    #includeFunc  flowRatePatch(name=OUTLET)"""

probe_templates["residual_probing"] = """
    #includeFunc residuals"""

probe_templates["inletoutletvelocity_probing"] = """
    #include "Probes_InletOutlet_Dict" """

probe_templates["xslice_probing"] = """
    #include "Probes_XSlices_Dict" """

probe_templates["midspanslice_probing"] = """
    #include "Probe_Slice_Dict" """

probe_templates["fieldave_probing"] = """
    #include "Probes_FieldAve_Dict" """

probe_templates["profile_probing"] = """
    #include "Probes_Profile_Dict" """

probe_templates["streamline_probing"] = """
    #include "Probes_Streamline_Dict" """


def get_template_contents():
    template_contents = {}
    for key, vallist in file_templates.items():
        template_contents[key] = {}
        for val in vallist:
            template_contents[key][val] = readtxtfile(os.path.join(key, val))

    return template_contents


