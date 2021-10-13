import os


def mgmt_parastud(settings):
    return 0


def mgmt_simulation(settings):
    print(settings)
    return 0


def create_jobmanagement(casetype, settings):
    case_templates = os.listdir(os.path.join(os.path.dirname(__file__), "../database/job_scriptsd"))

    if casetype == "parameterstudy":
        mgmt_parastud(settings)
    if casetype == "simulation":
        mgmt_simulation(settings)
    return 0
