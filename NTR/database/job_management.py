import os


import NTR


def mgmt_parastud(settings):
    return 0


def mgmt_simulation(settings):
    print(settings)
    return 0


def create_jobmanagement(casetype, settings):

    templatedir = os.path.join(os.path.dirname(NTR.__file__), "database", "job_scripts")
    case_templates = os.listdir(templatedir)
    job_settings = settings["case_settings"]["job"]

    scriptname = job_settings["job_script"]
    if casetype == "parameterstudy":
        mgmt_parastud(settings)
    if casetype == "simulation":
        mgmt_simulation(settings)
    return 0
