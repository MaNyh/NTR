import os
import shutil

import NTR
from NTR.database.case_dirstructure import casedirs
from NTR.utils.filehandling import write_file


def mgmt_parastud(settings, casepath):
    parastudsh_fname = "parastud.sh"
    txt = ""
    prepsh = "prep.sh"

    for idx,setting in enumerate(settings):
        subname = "paracase_" + str(idx)
        scriptname = setting["case_settings"]["job"]["job_script"]
        subcmd = setting["case_settings"]["sub_cmd"]
        txt += "# "+subname+ "\n"
        txt += "cd "+subname + "\n"
        txt += "sh " + prepsh+ "\n"
        txt += subcmd + " " + scriptname + "\n"
        txt += "cd ..\n"
    write_file(os.path.join(casepath, casedirs["simcase"], parastudsh_fname), txt)
    return 0


def mgmt_simulation(settings, scriptpath, casepath):
    targetpath = os.path.join(casepath, casedirs["simcase"], os.path.basename(scriptpath))
    shutil.copy(scriptpath, targetpath)

    job_procs = int(settings["case_settings"]["job"]["job_nodes"]) * int(settings["case_settings"]["job"]["job_ppn"])
    job_mempercpu = int(settings["case_settings"]["job"]["job_mem"]) / job_procs

    settings["case_settings"]["job"]["JOB_PROCS"] = str(job_procs)
    settings["case_settings"]["job"]["JOB_MEMPERCPU"] = str(job_mempercpu)
    for parameter in settings["case_settings"]["job"]:

        lookfor = parameter.upper()

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
    if "prep" not in settings["case_settings"].keys():
        return 0
    target = os.path.join(targetpath, prepfile)
    prep_text = settings["case_settings"]["prep"]
    write_file(target, prep_text)
    return 0


def write_runsim_bash(settings, targetpath):
    assert "sub_cmd" in settings["case_settings"].keys(), "sub_cmd (job-submission-command) not defined in settings"
    targetpath = os.path.join(targetpath, casedirs["simcase"])
    runfile = "runsim.sh"
    runcmd = settings["case_settings"]["sub_cmd"]
    text = "sh prep.sh"
    text += "\n" + runcmd + " " + settings["case_settings"]["job"]["job_script"] + ".sh"
    write_file(os.path.join(targetpath, runfile), text)


def create_jobmanagement(casetype, settings, casepath):
    templatedir = os.path.join(os.path.dirname(NTR.__file__), "database", "job_scripts")
    case_templates = os.listdir(templatedir)
    job_settings = settings["case_settings"]["job"]

    scriptname = job_settings["job_script"]
    scriptfile = scriptname + ".sh"
    assert scriptfile in case_templates, "chosen script is not available. allowed: " + str(case_templates)
    scriptpath = os.path.join(templatedir, scriptfile)

    # ToDo: this is stupid
    # if casetype == "parameterstudy":
    #    mgmt_parastud(settings, casepath)
    if casetype == "simulation":
        mgmt_simulation(settings, scriptpath, casepath)
    return 0
