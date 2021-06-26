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
                             "Probe_Slice_Dict",
                             "topoSetDict"],

                  "utils": ["monitor.py"],

                  ".": ["submit_job_hlrn.sh"],

                  }

probe_templates = {}

probe_templates["inletoutletfieldave_probing"] = """
    #include "Probes_inletoutletave_Dict" """

probe_templates["inletoutletvelocity_probing"] = """
    #include "Probes_InletOutlet_Dict" """

probe_templates["xslice_probing"] = """
    #include "Probes_XSlices_Dict" """

probe_templates["midspanslice_probing"] = """
    #include "Probe_Slice_Dict" """

probe_templates["profile_probing"] = """
    #include "Probes_Profile_Dict" """

probe_templates["streamline_probing"] = """
    #include "Probes_Streamline_Dict" """


