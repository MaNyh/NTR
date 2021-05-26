from NTR.preprocessing.openfoam.create_probes import createProbesProfileDict, createProbesStreamlineDict

#createProbesProfileDict("ressources/VTK/BLADE/BLADE_109000.vtk", 1, 1, 1, ".", tolerance=1e-6)
createProbesStreamlineDict("ressources/VTK/05_GWKVD_109000.vtk", 20, ".", 3, 135, 115, 0.0765)
