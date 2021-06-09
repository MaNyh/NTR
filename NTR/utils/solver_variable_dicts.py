openfoam_timeaveraged = {"p": "pMean",
                        "U": "UMean",
                        "T": "TMean",
                        "nut": "nutMean",
                        "alphat": "alphatMean",
                        "rho": "rhoMean"}

openfoam_instantanious = {"p": "p",
                            "U": "U",
                            "T": "T",
                            "nut": "nut",
                            "alphat": "alphat",
                            "rho": "rho"}

solver_var_dicts = {"openfoam_timeaveraged": openfoam_timeaveraged,
                    "openfoam_instantanious": openfoam_instantanious}
