import os
import shutil

import NTR
from NTR.database.case_dirstructure import casedirs
from NTR.utils.filehandling import write_file


def mgmt_parastud(settings, scriptname, casepath):
    return 0


def mgmt_simulation(settings, scriptpath, casepath):
    targetpath = os.path.join(casepath, casedirs["simcase"], os.path.basename(scriptpath))
    shutil.copy(scriptpath, targetpath)

    job_procs = int(settings["case_settings"]["job"]["job_nodes"]) * int(settings["case_settings"]["job"]["job_ppn"])
    job_mempercpu = int(settings["case_settings"]["job"]["job_mem"]) / job_procs

    settings["case_settings"]["job"]["JOB_PROCS"] = str(job_procs)
    settings["case_settings"]["job"]["JOB_MEMPERCPU"] = str(job_mempercpu)
    for parameter in settings["case_settings"]["job"]:
        if not parameter == "account":
            lookfor = parameter.upper()
        else:
            lookfor = parameter

        variable = settings["case_settings"]["job"][parameter]
        with open(targetpath) as fobj:
            newText = fobj.read().replace("<var " + lookfor + " var>", str(variable))
        with open(targetpath, "w") as fobj:
            fobj.write(newText)

    prep_target = os.path.join(casepath, casedirs["simcase"])
    write_prep_bash(settings, prep_target)

    return 0


def write_prep_bash(settings, targetpath):
    prepfile = "prep.sh"
    target = os.path.join(targetpath, prepfile)
    prep_text = settings["case_settings"]["prep"]
    write_file(target, prep_text)


def create_jobmanagement(casetype, settings, casepath):
    templatedir = os.path.join(os.path.dirname(NTR.__file__), "database", "job_scripts")
    case_templates = os.listdir(templatedir)
    job_settings = settings["case_settings"]["job"]

    scriptname = job_settings["job_script"]
    scriptfile = scriptname + ".sh"
    assert scriptfile in case_templates, "chosen script is not available"
    scriptpath = os.path.join(templatedir, scriptfile)

    if casetype == "parameterstudy":
        mgmt_parastud(settings, scriptpath, casepath)
    if casetype == "simulation":
        mgmt_simulation(settings, scriptpath, casepath)
    return 0
