from NTR.preprocessing.openfoam.create_probes import createProbesProfileDict, createProbesStreamlineDict

createProbesProfileDict("ressources/VTK/BLADE/BLADE_600.vtk", 1, 1, 1, ".", tolerance=1e-6)
