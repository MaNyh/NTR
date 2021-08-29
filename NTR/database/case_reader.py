import os
import re
from functools import reduce  # forward compatibility for Python 3
import operator

from NTR.utils.filehandling import get_directory_structure, write_yaml_dict, yaml_dict_read ,write_pickle, read_pickle

def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)

def setInDict(dataDict, mapList, value):
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value

def nested_dict_pairs_iterator(dict_obj):
    ''' This function accepts a nested dictionary as argument
        and iterate over all values of nested dictionaries
    '''
    # Iterate over all key-value pairs of dict argument
    for key, value in dict_obj.items():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over all its values
            for pair in  nested_dict_pairs_iterator(value):
                yield (key, *pair)
        else:
            # If value is not dict type then yield the value
            yield (key, value)



def find_vars_opts(case_structure):
    # allowing names like JOB_NUMBERS, only capital letters and underlines - no digits, no whitespaces
    varsignature = r"<var[\s[A-Z]{3,}_{,}[A-Z]{3,}]{,}\svar>"
    optsignature = r"<opt[\s[A-Z]{3,}_{,}[A-Z]{3,}]{,}\sopt>"
    siglim = (5, -5)

    parameters = []
    options = []
    all_pairs = list(nested_dict_pairs_iterator(case_structure))
    for pair in all_pairs:
        filepath = os.path.join(*pair[:-1])
        with open(os.path.join(os.path.dirname(__file__), "case_templates", filepath), "r") as fhandle:
            for line in fhandle.readlines():

                lookup_var = re.search(varsignature, line)
                if lookup_var:
                    span = lookup_var.span()
                    parameter = line[span[0] + siglim[0]:span[1] + siglim[1]]
                    setInDict(case_structure,pair[:-1],{parameter:"var"})
                lookup_opt = re.search(optsignature, line)
                if lookup_opt:
                    span = lookup_opt.span()
                    opt_name = line[span[0] + siglim[0]:span[1] + siglim[1]]
                    setInDict(case_structure,pair[:-1],{opt_name:"opt"})

    return case_structure


def create_simulationcase(path_to_yaml_dict):
    casedirectories = {"ressources": "00_Ressources",
                       "meshing": "01_Meshing",
                       "simcase": "02_Simcase",
                       "solution": "03_Solution",
                       "data": "04_Data"}

    case_templates = os.listdir(os.path.join(os.path.dirname(__file__), "case_templates"))

    case_structures = {}
    for case_name in case_templates:
        cstruct = get_directory_structure(
            os.path.join(os.path.dirname(__file__), "case_templates", case_name))
        case_structures[case_name] = cstruct

    settings_dict = yaml_dict_read(path_to_yaml_dict)

    assert "name" in settings_dict["case_settings"], "no name for the case defined"
    casename = settings_dict["case_settings"]["name"]

    settings = yaml_dict_read(path_to_yaml_dict)
    casepath = os.path.abspath(os.path.dirname(path_to_yaml_dict))

    case_type = settings["case_settings"]["case_type"]
    case_directories = case_structures[case_type]
    assert case_type in case_structures.keys(), "case_type " + case_type + " not found in templates."

    path_to_geo_ressources = os.path.join(casepath, "04_Data", "geometry.pkl")
    assert os.path.isfile(path_to_geo_ressources), "no geometry.pkl found, create the geometry first"
    geo_ressources = read_pickle(os.path.join(path_to_geo_ressources))

    create_casedirstructure(casedirectories,casepath)
    case_structure = find_vars_opts(case_structures[case_type])
    necessarities = list_necessary_parameters(case_structure)


def list_necessary_parameters(case_structure):

    return []

def get_parametrized_simstructure(template_path):
    case_structures = {}
    case_structures[case_name] = get_directory_structure(template_path)


def create_casedirstructure(casedirectories,casepath):
    for d in casedirectories.values():
        if not os.path.isdir(os.path.join(casepath, d)):
            os.mkdir(os.path.join(casepath, d))

def create_simdirstructure(filetemplates,path):
    directories = list(filetemplates.keys())
    for d in directories:
        if not os.path.isdir(os.path.join(path, d)):
            os.mkdir(os.path.join(path, d))
    return 0

def test_create_simulationcase(tmpdir):
    test_dict = {"case_settings": {"case_type": "openfoam_cascade_les","name":"testcase"}}
    test_file = tmpdir / "test_create_simulationcase.yml"
    test_geo_dict = {}
    os.mkdir(tmpdir/"04_Data")
    test_geo_file = tmpdir / os.path.join("04_Data", "geometry.pkl")
    write_yaml_dict(test_file,test_dict)
    write_pickle(test_geo_file,test_geo_dict)
    create_simulationcase(test_file)
    print(tmpdir)
