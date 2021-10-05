# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import os
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CGNS Series Reader'
tRACEcgns = CGNSSeriesReader(FileNames=["tmp.cgns"])

# create a new 'Merge Blocks'
mergeBlocks1 = MergeBlocks(Input=tRACEcgns)

# save data
SaveData('solution.vtk', proxy=mergeBlocks1, FileType='Ascii')

