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

# Properties modified on tRACEcgns
tRACEcgns.CellArrayStatus = ['Density', 'Mach', 'Pressure', 'Temperature', 'TurbulentDissipationRate', 'TurbulentDistance', 'TurbulentEnergyKinetic', 'Velocity', 'ViscosityEddy', 'ViscosityEddyRatio', 'diffusionCoeff', 'dissipationCoeff', 'productionCoeff']
# save data
SaveData(r"<var BLADE var>", proxy=tRACEcgns)
"""
# Properties modified on tRACEcgns
tRACEcgns.Blocks = ['/Hierarchy/Base#1/Main_Blade_down.1/3', '/Hierarchy/Base#1/Main_Blade_down.2/5', '/Hierarchy/Base#1/Main_Blade_hubgap1.1/3', '/Hierarchy/Base#1/Main_Blade_hubgap1.2/5', '/Hierarchy/Base#1/Main_Blade_hubgap1.3/7', '/Hierarchy/Base#1/Main_Blade_hubgap2/3', '/Hierarchy/Base#1/Main_Blade_inlet/3', '/Hierarchy/Base#1/Main_Blade_outlet/3', '/Hierarchy/Base#1/Main_Blade_skin.1/3', '/Hierarchy/Base#1/Main_Blade_skin.2/4', '/Hierarchy/Base#1/Main_Blade_skin.3/5', '/Hierarchy/Base#1/Main_Blade_up.1/3', '/Hierarchy/Base#1/Main_Blade_up.2/5', '/Hierarchy/Base#1/domain1.01/3', '/Hierarchy/Base#1/domain1.03/4', '/Hierarchy/Base#1/domain1.04/4', '/Hierarchy/Base#1/domain1.05/5', '/Hierarchy/Base#1/domain1.06/5', '/Hierarchy/Base#1/domain1.07/6', '/Hierarchy/Base#1/domain3.01/3', '/Hierarchy/Base#1/domain3.03/4', '/Hierarchy/Base#1/domain3.04/4', '/Hierarchy/Base#1/domain3.05/5']
# save data
SaveData(r"<var HUB var>", proxy=tRACEcgns)
"""

# Properties modified on tRACEcgns
tRACEcgns.Blocks = ['/Hierarchy/Base#1/domain1.01/5', '/Hierarchy/Base#1/domain1.02/6', '/Hierarchy/Base#1/domain1.04/6', '/Hierarchy/Base#1/domain1.08/7', '/Hierarchy/Base#1/domain1.10/7', '/Hierarchy/Base#1/domain1.14/8']
# save data
SaveData(r"<var OUTLET var>", proxy=tRACEcgns)
# Properties modified on tRACEcgns

"""
tRACEcgns.Blocks = ['/Hierarchy/Base#1/Main_Blade_down.3/4', '/Hierarchy/Base#1/Main_Blade_down.4/6', '/Hierarchy/Base#1/Main_Blade_inlet/4', '/Hierarchy/Base#1/Main_Blade_outlet/4', '/Hierarchy/Base#1/Main_Blade_skin1.3/4', '/Hierarchy/Base#1/Main_Blade_skin1.5/5', '/Hierarchy/Base#1/Main_Blade_skin1.7/6', '/Hierarchy/Base#1/Main_Blade_skin1.8/8', '/Hierarchy/Base#1/Main_Blade_skin1.9/7', '/Hierarchy/Base#1/Main_Blade_up.3/4', '/Hierarchy/Base#1/Main_Blade_up.4/6', '/Hierarchy/Base#1/domain1.08/4', '/Hierarchy/Base#1/domain1.13/5', '/Hierarchy/Base#1/domain1.14/5', '/Hierarchy/Base#1/domain1.15/6', '/Hierarchy/Base#1/domain1.16/6', '/Hierarchy/Base#1/domain3.06/4', '/Hierarchy/Base#1/domain3.10/5', '/Hierarchy/Base#1/domain3.11/5']
# save data
SaveData(r"<var SHROUD var>", proxy=tRACEcgns)

"""
# Properties modified on tRACEcgns
tRACEcgns.Blocks = ['/Hierarchy/Base#1/domain3.04/6', '/Hierarchy/Base#1/domain3.05/7', '/Hierarchy/Base#1/domain3.08/7', '/Hierarchy/Base#1/domain3.09/8', '/Hierarchy/Base#1/domain3.10/9', '/Hierarchy/Base#1/domain3.11/8']
# save data
SaveData(r"<var INLET var>", proxy=tRACEcgns)


# Properties modified on tRACEcgns
tRACEcgns.Blocks = ['/Hierarchy/Base#1/Main_Blade_hubgap1.1/4', '/Hierarchy/Base#1/Main_Blade_hubgap1.2/6', '/Hierarchy/Base#1/Main_Blade_hubgap1.3/8', '/Hierarchy/Base#1/Main_Blade_hubgap2/4']
# save data
SaveData(r"<var TIP var>", proxy=tRACEcgns)
