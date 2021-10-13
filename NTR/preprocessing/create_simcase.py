import os
import re
import shutil
import itertools
import copy
import tempfile
import yaml
import glob

from NTR.utils.dicthandling import setInDict, nested_val_set, nested_dict_pairs_iterator
from NTR.utils.filehandling import get_directory_structure, yaml_dict_read, read_pickle
from NTR.utils.functions import func_by_name
from NTR.database.case_dirstructure import casedirs


def find_vars_opts(case_structure):
    # allowing names like JOB_NUMBERS, only capital letters and underlines - no digits, no whitespaces
    varsignature = r"<var [A-Z]{3,}(_{1,1}[A-Z]{3,}){,} var>"
    optsignature = r"<opt [A-Z]{3,}(_{1,1}[A-Z]{3,}){,} opt>"
    siglim = (5, -5)

    all_pairs = list(nested_dict_pairs_iterator(case_structure))
    for pair in all_pairs:
        setInDict(case_structure, pair[:-1], {})
        filepath = os.path.join(*pair[:-1])
        with open(os.path.join(os.path.dirname(__file__), "../database/case_templates", filepath), "r") as fhandle:
            for line in fhandle.readlines():

                lookforvar = True
                lookforopt = True

                while (lookforvar):
                    lookup_var = re.search(varsignature, line)
                    if not lookup_var:
                        lookforvar = False
                    else:
                        span = lookup_var.span()
                        parameter = line[span[0] + siglim[0]:span[1] + siglim[1]]
                        setInDict(case_structure, list(pair[:-1]) + [parameter], "var")
                        match = line[span[0]:span[1]]
                        line = line.replace(match, "")

                while (lookforopt):
                    lookup_opt = re.search(optsignature, line)
                    if not lookup_opt:
                        lookforopt = False
                    else:
                        span = lookup_opt.span()
                        opt_name = line[span[0] + siglim[0]:span[1] + siglim[1]]
                        setInDict(case_structure, list(pair[:-1]) + [opt_name], "opt")
                        match = line[span[0]:span[1]]
                        line = line.replace(match, "")
    return case_structure


def read_parastudyaml(path_to_yaml_dict):
    settings = yaml_dict_read(path_to_yaml_dict)

    allvals = list(nested_dict_pairs_iterator(settings))
    allvals_aslist = [[i[-1]] if type(i[-1]) != list else i[-1] for i in allvals]
    allstruct = [i[:-1] for i in allvals]

    sets = list(itertools.product(*allvals_aslist))

    settings_parastud = []

    for vals in sets:
        kwargs = settings.copy()
        for idx, dict_struct in enumerate(allstruct):
            val = vals[idx]
            nested_val_set(kwargs, dict_struct, val)
        settings_parastud.append(copy.deepcopy(kwargs))

    return settings_parastud


def mgmt_parastud(settings):
    return 0


def mgmt_simulation(settings):
    return 0


def create_jobmanagement(casetype, settings):
    if casetype == "parameterstudy":
        mgmt_parastud(settings)
    if casetype == "simulation":
        mgmt_simulation(settings)
    return 0


def create_parastudsims(path_to_parayaml):
    yamldict = yaml_dict_read(path_to_parayaml)
    casetype = yamldict["case_settings"]["type"]
    assert casetype == "parameterstudy", "check your yaml-dict. the case is not defined as a parameterstudy"

    settings = read_parastudyaml(path_to_parayaml)

    casepath = os.path.abspath(os.path.dirname(path_to_parayaml))
    sim_dirs = []
    for idx, settings_dict in enumerate(settings):
        settings_dict["case_settings"]["type"] = "simulation"
        subname = "paracase_" + str(idx)
        tmp_dir = tempfile.TemporaryDirectory()
        target_dir = os.path.join(casepath, subname)
        tmp_yml = os.path.join(tmp_dir.name, "tmp_settings.yaml")
        with open(tmp_yml, "w") as handle:
            yaml.dump(settings_dict, handle, default_flow_style=False)
        create_simulationcase(tmp_yml, subname)
        files = glob.glob(os.path.join(tmp_dir.name, "02_Simcase") + "/*")
        for f in files:
            shutil.move(f, target_dir)
        shutil.copyfile(tmp_yml, os.path.join(target_dir, subname + "_settings.yml"))
        tmp_dir.cleanup()
        sim_dirs.append(target_dir)

    create_jobmanagement(casetype, settings)


def create_simulationcase(path_to_yaml_dict, subdir=False):
    case_templates = os.listdir(os.path.join(os.path.dirname(__file__), "../database/case_templates"))

    case_structures = {}
    for cname in case_templates:
        cstruct = get_directory_structure(
            os.path.join(os.path.dirname(__file__), "../database/case_templates", cname))
        case_structures[cname] = cstruct

    settings = yaml_dict_read(path_to_yaml_dict)
    casetype = yamldict["case_settings"]["type"]
    assert casetype == "simulation", "check your yaml-dict. the case is not defined as a simulation"

    assert "name" in settings["case_settings"], "no name for the case defined"

    casepath = os.path.abspath(os.path.dirname(path_to_yaml_dict))

    case_type = settings["case_settings"]["case_type"]
    assert case_type in case_structures.keys(), "case_type " + case_type + " not found in templates."

    if subdir == False:
        path_to_sim = os.path.join(casepath, casedirs["simcase"])
    else:
        path_to_sim = os.path.join(casepath, casedirs["simcase"], subdir)
        os.mkdir(os.path.join(casepath, casedirs["simcase"]))
        os.mkdir(path_to_sim)

    create_casedirstructure(casedirs, casepath)
    case_structure = case_structures[case_type]
    create_simdirstructure(case_structure, path_to_sim)
    copy_template(case_type, case_structure, path_to_sim)
    case_structure_parameters = find_vars_opts(case_structure)
    check_settings_necessarities(case_structure_parameters, settings)
    writeout_simulation(case_structure_parameters, path_to_sim, settings)
    writeout_simulation_options(case_structure_parameters, path_to_sim, settings)


def writeout_simulation(case_structure_parameters, path_to_sim, settings):
    walk_casefile_list = nested_dict_pairs_iterator(case_structure_parameters)
    for parameterdata in walk_casefile_list:
        fpath = os.path.join(path_to_sim, *parameterdata[1:-2])
        parametername = parameterdata[-2]

        para_type = parameterdata[-1]
        if para_type == "var":
            para_type = "variables"
            variable = settings["simcase_settings"][para_type][parametername]
            with open(fpath) as fobj:
                newText = fobj.read().replace("<var " + parametername + " var>", str(variable))
            with open(fpath, "w") as fobj:
                fobj.write(newText)


def writeout_simulation_options(case_structure_parameters, path_to_sim, settings):
    walk_casefile_list = nested_dict_pairs_iterator(case_structure_parameters)
    for parameterdata in walk_casefile_list:
        fpath = os.path.join(path_to_sim, *parameterdata[1:-2])
        parametername = parameterdata[-2]

        para_type = parameterdata[-1]
        if para_type == "opt":
            para_type = "options"
            option = settings["simcase_settings"][para_type][parametername]
            if option == True:
                assert "simcase_optiondef" in list(settings.keys()), "simcase_optiondef not defined in configuration"
                assert parametername in list(
                    settings["simcase_optiondef"].keys()), parametername + " not defined in simcase_optiondef"
                optiondefinition = settings["simcase_optiondef"][parametername]

                optionfunc = func_by_name(optiondefinition["func"])
                optionargs = optiondefinition["args"]
                optionargs["path_to_sim"] = path_to_sim
                optionargs["case_settings"] = settings
                geo_fpath = os.path.join(path_to_sim, "..", "04_Data", "geometry.pkl")
                if os.path.isfile(geo_fpath):
                    optionargs["geomdat_dict"] = read_pickle(os.path.join(path_to_sim, "..", "04_Data", "geometry.pkl"))
                else:
                    optionargs["geomdat_dict"] = None
                    print("warning, no geometry.pkl found. uncertain if this is going to work...")
                optionfunc(**optionargs)

                replace_opt = optiondefinition["insert"]
                with open(fpath) as fobj:
                    newText = fobj.read().replace("<opt " + parametername + " opt>", str(replace_opt))
                with open(fpath, "w") as fobj:
                    fobj.write(newText)
            else:
                with open(fpath) as fobj:
                    newText = fobj.read().replace("<opt " + parametername + " opt>", "")
                with open(fpath, "w") as fobj:
                    fobj.write(newText)


def copy_template(case_type, case_structure, path_to_sim):
    for file in nested_dict_pairs_iterator(case_structure):
        filename = file[-2]
        dirstructure = file[1:-2]
        if dirstructure == ():
            dirstructure = ""

        template_fpath = os.path.join(os.path.dirname(__file__), "../database/case_templates", case_type, *dirstructure,
                                      filename)
        sim_fpath = os.path.join(path_to_sim, *dirstructure, filename)
        shutil.copyfile(template_fpath, sim_fpath)


def check_settings_necessarities(case_structure, settings_dict):
    necessarities = list(nested_dict_pairs_iterator(case_structure))
    necessarity_vars = []
    necessarity_opts = []
    for item in necessarities:
        if item[-1] == "var":
            necessarity_vars.append(item[-2])
        if item[-1] == "opt":
            necessarity_opts.append(item[-2])

    assert "variables" in settings_dict["simcase_settings"].keys(), "variables not set simcase_settings"
    if settings_dict["simcase_settings"]["variables"]:
        settings_variables = list(settings_dict["simcase_settings"]["variables"].keys())
    else:
        settings_variables = []

    assert "options" in settings_dict["simcase_settings"].keys(), "options not set in simcase_settings"
    if settings_dict["simcase_settings"]["options"]:
        settings_options = list(settings_dict["simcase_settings"]["options"].keys())
    else:
        settings_options = []

    for variable in necessarity_vars:
        assert variable in settings_variables, "variable " + variable + " not set in configuration file"
    for option in necessarity_opts:
        assert option in settings_options, "option " + option + " not set in configuration file"


def get_parametrized_simstructure(template_path):
    case_structures = {}
    case_structures[case_name] = get_directory_structure(template_path)


def create_casedirstructure(casedirectories, casepath):
    for d in casedirectories.values():
        if not os.path.isdir(os.path.join(casepath, d)):
            os.mkdir(os.path.join(casepath, d))


def create_simdirstructure(filetemplates, path):
    directories = list(nested_dict_pairs_iterator(filetemplates))
    for d in directories:
        dirstructure = d[1:-2]
        if dirstructure == ():
            dirstructure = ""
        for dir in dirstructure:
            dpath = os.path.join(path, dir)
            if not os.path.isdir(dpath):
                os.mkdir(dpath)
    return 0
