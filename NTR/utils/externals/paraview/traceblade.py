# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CGNS Series Reader'
tRACE_2Dcgns = CGNSSeriesReader(FileNames=[r"<var MESH var>"])
# Properties modified on tRACE_2Dcgns
tRACE_2Dcgns.CellArrayStatus = ['BoundaryLayerEdgeIndex', 'BoundaryLayerEdgeVelocity', 'Density', 'DistanceWallCoordinate', 'HeatFluxWall', 'HeightLaminarSublayer', 'Intermittency', 'Pressure', 'ReynoldsMomentumThickness', 'ReynoldsMomentumThickness0', 'ShearStressWall', 'Temperature', 'ThicknessBoundaryLayer', 'ThicknessDisplacement', 'ThicknessDisplacementInc', 'ThicknessMomentum', 'ThicknessMomentumInc', 'TurbulenceIntensityNearWall', 'TurbulentDissipationRate', 'TurbulentEnergyKinetic', 'Velocity', 'ViscosityEddyRatioNearWall']



# Properties modified on tRACE_2Dcgns
tRACE_2Dcgns.Blocks = ['/Hierarchy/Base#1/OUTLET-z08-p0/Grid', '/Hierarchy/Base#1/OUTLET-z09-p0/Grid', '/Hierarchy/Base#1/OUTLET-z10-p0/Grid']


# create a new 'Merge Blocks'
mergeBlocks1 = MergeBlocks(Input=tRACE_2Dcgns)


# save data
SaveData(r"<var OUTLET var>", proxy=mergeBlocks1)

# Properties modified on tRACE_2Dcgns
tRACE_2Dcgns.Blocks = ['/Hierarchy/Base#1/INLET-z01-p0/Grid', '/Hierarchy/Base#1/INLET-z05-p0/Grid', '/Hierarchy/Base#1/INLET-z06-p0/Grid']
# save data
SaveData(r"<var INLET var>", proxy=mergeBlocks1)

# Properties modified on tRACE_2Dcgns
tRACE_2Dcgns.Blocks = ['/Hierarchy/Base#1/BLADE-z02-p0/Grid', '/Hierarchy/Base#1/BLADE-z03-p0/Grid', '/Hierarchy/Base#1/BLADE-z04-p0/Grid', '/Hierarchy/Base#1/BLADE-z12-p0/Grid']
# save data
SaveData(r"<var BLADE var>", proxy=mergeBlocks1)
