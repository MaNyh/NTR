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
tRACE_2Dcgns.Blocks = ['/Hierarchy/Base#1/outlet-z27-p1/Grid', '/Hierarchy/Base#1/outlet-z28-p0/Grid', '/Hierarchy/Base#1/outlet-z30-p1/Grid', '/Hierarchy/Base#1/outlet-z34-p1/Grid', '/Hierarchy/Base#1/outlet-z36-p0/Grid', '/Hierarchy/Base#1/outlet-z40-p1/Grid']

# create a new 'Merge Blocks'
mergeBlocks1 = MergeBlocks(Input=tRACE_2Dcgns)


# save data
SaveData(r"<var OUTLET var>", proxy=mergeBlocks1)

# Properties modified on tRACE_2Dcgns
tRACE_2Dcgns.Blocks = ['/Hierarchy/Base#1/inlet-z46-p1/Grid', '/Hierarchy/Base#1/inlet-z47-p1/Grid', '/Hierarchy/Base#1/inlet-z50-p0/Grid', '/Hierarchy/Base#1/inlet-z51-p0/Grid', '/Hierarchy/Base#1/inlet-z52-p1/Grid', '/Hierarchy/Base#1/inlet-z53-p1/Grid']

# save data
SaveData(r"<var INLET var>", proxy=mergeBlocks1)

# Properties modified on tRACE_2Dcgns
tRACE_2Dcgns.Blocks = ['/Hierarchy/Base#1/blade-z14-p0/Grid', '/Hierarchy/Base#1/blade-z15-p0/Grid', '/Hierarchy/Base#1/blade-z16-p0/Grid', '/Hierarchy/Base#1/blade-z17-p0/Grid', '/Hierarchy/Base#1/blade-z18-p0/Grid', '/Hierarchy/Base#1/blade-z19-p0/Grid', '/Hierarchy/Base#1/blade-z20-p0/Grid', '/Hierarchy/Base#1/blade-z21-p0/Grid', '/Hierarchy/Base#1/blade-z22-p0/Grid']
# save data
SaveData(r"<var BLADE var>", proxy=mergeBlocks1)
