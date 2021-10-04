# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:33:27 2019

@author: Mark Ziesse / Malte Nyhuis
"""

"""
def openFoam_create_of_les_probe_dicts(case_settings, geo_ressources):

    midspan_z = case_settings["mesh"]["extrudeLength"] / 2

    output_path = os.path.join("02_Simcase", "system")

    probes = {}

    if case_settings["probing"]["probes"]["profile_probing"]:
        sampling_rate = case_settings["probes"]["profile_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        outprobes = openFoam_createProbesProfileDict(geo_ressources,
                                                     midspan_z,
                                                     case_settings["probes"]["profile_probes"]["pden_ps"],
                                                     case_settings["probes"]["profile_probes"]["pden_ss"],
                                                     timestepinterval,
                                                     output_path,
                                                     case_settings["probes"]["profile_probes"]["start_time"],
                                                     case_settings["probes"]["profile_probes"]["end_time"],
                                                     case_settings["probes"]["profile_probes"]["tolerance"]
                                                     )
        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["streamline_probing"]:

        sampling_rate = case_settings["probes"]["streamline_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        outprobes = openFoam_createProbesStreamlineDict(case_settings["probes"]["streamline_probes"]["nop_streamline"],
                                                        output_path, timestepinterval,
                                                        case_settings["probes"]["streamline_probes"]["start_time"],
                                                        case_settings["probes"]["streamline_probes"]["end_time"],
                                                        geo_ressources)

        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["inletoutletvelocity_probing"]:
        sampling_rate = case_settings["probes"]["inletoutlet_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        outprobes = openFoam_createProbesInletOutlet(geo_ressources, timestepinterval, output_path,
                                                     case_settings["probes"]["inletoutlet_probes"]["start_time"],
                                                     case_settings["probes"]["inletoutlet_probes"]["end_time"])
        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["xslice_probing"]:

        sampling_rate = case_settings["probes"]["xsclicing_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        outprobes = openFoam_createXSliceProbes(geo_ressources, case_settings["probes"]["xsclicing_probes"]["nop"],
                                                case_settings["probes"]["xsclicing_probes"]["x_slice_one"],
                                                case_settings["probes"]["xsclicing_probes"]["x_slice_two"],
                                                timestepinterval, output_path,
                                                case_settings["probes"]["xsclicing_probes"]["start_time"],
                                                case_settings["probes"]["xsclicing_probes"]["end_time"])
        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["inletoutletfieldave_probing"]:
        sampling_rate = case_settings["probes"]["inletoutletfieldave_probing"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        openFoam_create_inletoutletave_probe_dict(case_settings["probes"]["inletoutletfieldave_probing"]["start_time"],
                                                  case_settings["probes"]["inletoutletfieldave_probing"]["end_time"],
                                                  timestepinterval,
                                                  output_path)

    if case_settings["probing"]["probes"]["vk_stagnationflow_probing"]:
        sampling_rate = case_settings["probes"]["inletoutletfieldave_probing"]["sampling_rate"]

        U = [float(i) for i in case_settings["case"]["case_parameters"]["U"]["UINLET"].split(" ")]
        angle = np.arccos(U[0] / U[1]) * 180 / np.pi
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))

        outprobes = openFoam_create_vk_stagflow_probes(geo_ressources,
                                                       case_settings["probes"]["vk_stagnationflow_probing"]["nop"],
                                                       case_settings["probes"]["vk_stagnationflow_probing"]["length"],
                                                       angle, timestepinterval, output_path,
                                                       case_settings["probes"]["inletoutletfieldave_probing"][
                                                           "start_time"],
                                                       case_settings["probes"]["inletoutletfieldave_probing"][
                                                           "end_time"])
        for k, v in outprobes.items():
            probes[k] = v

    #midspan_z = geo_ressources["span_z"]/2
    y_lower = geo_ressources["periodics"]["ylower"]
    y_upper = geo_ressources["periodics"]["yupper"]
    ssPoly = geo_ressources["sidePolys"]["ssPoly"]
    psPoly = geo_ressources["sidePolys"]["psPoly"]
    inlet = geo_ressources["flowbounds"]["inletPoly"]
    outlet = geo_ressources["flowbounds"]["outletPoly"]

    plt.close('all')
    plt.figure(figsize=(8, 8))
    plt.plot(inlet.points[::,0],inlet.points[::,1], '-r', lw=1, label='inlet')
    plt.plot(outlet.points[::,0],outlet.points[::,1], '-b', lw=1, label='outlet')
    plt.plot(ssPoly.points[::,0],ssPoly.points[::,1], '-k', lw=1, label='ss_profil')
    plt.plot(psPoly.points[::,0],psPoly.points[::,1], '-k', lw=1, label='ps_profil')
    plt.plot(y_lower.points[::,0],y_lower.points[::,1], '-y', lw=1, label='lower_peri')
    plt.plot(y_upper.points[::,0],y_upper.points[::,1], '-c', lw=1, label='upper_peri')

    allowed_colors = list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())
    del allowed_colors[allowed_colors.index("white")]
    del allowed_colors[allowed_colors.index("whitesmoke")]
    del allowed_colors[allowed_colors.index("seashell")]
    del allowed_colors[allowed_colors.index("linen")]
    del allowed_colors[allowed_colors.index("floralwhite")]
    del allowed_colors[allowed_colors.index("azure")]
    del allowed_colors[allowed_colors.index("mintcream")]
    del allowed_colors[allowed_colors.index("ghostwhite")]
    del allowed_colors[allowed_colors.index("aliceblue")]
    del allowed_colors[allowed_colors.index("lavenderblush")]
    del allowed_colors[allowed_colors.index("lightcyan")]
    del allowed_colors[allowed_colors.index("oldlace")]
    del allowed_colors[allowed_colors.index("antiquewhite")]
    del allowed_colors[allowed_colors.index("mistyrose")]

    for k, v in probes.items():
        color = allowed_colors.pop(np.random.randint(len(allowed_colors) - 1))

        plt.plot(v[0], v[1], color, marker="X", lw=1, markersize=3, label=k)

    plt.legend(loc='best')
    plt.savefig(os.path.join('kontrollplot_probes.pdf'))
    plt.close('all')
"""

