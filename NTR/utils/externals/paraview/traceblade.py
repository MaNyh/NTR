# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CGNS Series Reader'
tRACEcgns = CGNSSeriesReader(FileNames=[r"<var MESH var>"])
tRACEcgns.Blocks = ['/Hierarchy/Base#1/Block_1/Grid', '/Hierarchy/Base#1/Block_10/Grid', '/Hierarchy/Base#1/Block_11/Grid', '/Hierarchy/Base#1/Block_12/Grid', '/Hierarchy/Base#1/Block_2/Grid', '/Hierarchy/Base#1/Block_3/Grid', '/Hierarchy/Base#1/Block_4/Grid', '/Hierarchy/Base#1/Block_5/Grid', '/Hierarchy/Base#1/Block_6/Grid', '/Hierarchy/Base#1/Block_7/Grid', '/Hierarchy/Base#1/Block_8/Grid', '/Hierarchy/Base#1/Block_9/Grid']

# Properties modified on tRACEcgns
tRACEcgns.CellArrayStatus = ['Density', 'Mach', 'Pressure', 'Temperature', 'TurbulentDissipationRate', 'TurbulentDistance', 'TurbulentEnergyKinetic', 'Velocity', 'ViscosityEddy', 'ViscosityEddyRatio', 'diffusionCoeff', 'dissipationCoeff', 'productionCoeff']

# Properties modified on tRACEcgns
tRACEcgns.Blocks = ['/Hierarchy/Base#1/Block_5/2', '/Hierarchy/Base#1/Block_6/2', '/Hierarchy/Base#1/Block_7/2']
# save data
SaveData(r"<var OUTLET var>", proxy=tRACEcgns)

# Properties modified on tRACEcgns
tRACEcgns.Blocks = ['/Hierarchy/Base#1/Block_1/1', '/Hierarchy/Base#1/Block_2/1', '/Hierarchy/Base#1/Block_3/1']
# save data
SaveData(r"<var INLET var>", proxy=tRACEcgns)

# Properties modified on tRACEcgns
tRACEcgns.Blocks = ['/Hierarchy/Base#1/Block_10/4', '/Hierarchy/Base#1/Block_11/4', '/Hierarchy/Base#1/Block_12/3', '/Hierarchy/Base#1/Block_9/2']
# save data
SaveData(r"<var BLADE var>", proxy=tRACEcgns)
