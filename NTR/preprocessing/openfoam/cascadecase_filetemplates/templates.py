import os


def readtxtfile(path_to_file):
    path = os.path.abspath(os.path.dirname(__file__))
    fpath = os.path.join(path, path_to_file)
    with open(fpath, "r") as fobj:
        content = fobj.readlines()
    return "".join(content)


file_templates = {"0": ["alphat", "nut", "p", "T", "U"],
                  "constant": ["thermophysicalProperties", "turbulenceProperties"],
                  "system": ["controlDict", "createPatchDict", "decomposeParDict", "fvSchemes", "fvSolution",
                             "mapFieldsDict", "SurfaceSampleDict", "topoSetDict"],
                  "utils": ["monitor.py"]
                  }


def get_template_contents():
    template_contents = {}
    for key, vallist in file_templates.items():
        template_contents[key] = {}
        for val in vallist:
            template_contents[key][val] = readtxtfile(os.path.join("cascadecase_filetemplates", key, val))

    return template_contents
