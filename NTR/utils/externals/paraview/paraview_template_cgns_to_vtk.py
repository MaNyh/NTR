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
#tRACEcgns.Blocks = ['/Hierarchy/Base#1/Block_1/Grid', '/Hierarchy/Base#1/Block_10/Grid', '/Hierarchy/Base#1/Block_11/Grid', '/Hierarchy/Base#1/Block_12/Grid', '/Hierarchy/Base#1/Block_2/Grid', '/Hierarchy/Base#1/Block_3/Grid', '/Hierarchy/Base#1/Block_4/Grid', '/Hierarchy/Base#1/Block_5/Grid', '/Hierarchy/Base#1/Block_6/Grid', '/Hierarchy/Base#1/Block_7/Grid', '/Hierarchy/Base#1/Block_8/Grid', '/Hierarchy/Base#1/Block_9/Grid']
#tRACEcgns.CellArrayStatus = []

# get active view
#renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1436, 616]

# show data in view
#tRACEcgnsDisplay = Show(tRACEcgns, renderView1)

# trace defaults for the display properties.
#tRACEcgnsDisplay.Representation = 'Outline'
#tRACEcgnsDisplay.ColorArrayName = ['POINTS', '']
#tRACEcgnsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#tRACEcgnsDisplay.SelectOrientationVectors = 'None'
#tRACEcgnsDisplay.ScaleFactor = 0.03700000000000002
#tRACEcgnsDisplay.SelectScaleArray = 'None'
#tRACEcgnsDisplay.GlyphType = 'Arrow'
#tRACEcgnsDisplay.GlyphTableIndexArray = 'None'
#tRACEcgnsDisplay.GaussianRadius = 0.001850000000000001
#tRACEcgnsDisplay.SetScaleArray = [None, '']
#tRACEcgnsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#tRACEcgnsDisplay.OpacityArray = [None, '']
#tRACEcgnsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
#tRACEcgnsDisplay.DataAxesGrid = 'GridAxesRepresentation'
#tRACEcgnsDisplay.PolarAxes = 'PolarAxesRepresentation'
#tRACEcgnsDisplay.ScalarOpacityUnitDistance = 0.0065557524973083444

# reset view to fit data
#renderView1.ResetCamera()

# get the material library
#materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
#renderView1.Update()

# set scalar coloring
#ColorBy(tRACEcgnsDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
#tRACEcgnsDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
#vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
#vtkBlockColorsLUT.InterpretValuesAsCategories = 1
#vtkBlockColorsLUT.AnnotationsInitialized = 1
#vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
#vtkBlockColorsLUT.ActiveAnnotatedValues = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
#vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# get opacity transfer function/opacity map for 'vtkBlockColors'
#vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# create a new 'Merge Blocks'
mergeBlocks1 = MergeBlocks(Input=tRACEcgns)

# show data in view
#mergeBlocks1Display = Show(mergeBlocks1, renderView1)

# trace defaults for the display properties.
#mergeBlocks1Display.Representation = 'Surface'
#mergeBlocks1Display.ColorArrayName = [None, '']
#mergeBlocks1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#mergeBlocks1Display.SelectOrientationVectors = 'None'
#mergeBlocks1Display.ScaleFactor = 0.03700000000000002
#mergeBlocks1Display.SelectScaleArray = 'None'
#mergeBlocks1Display.GlyphType = 'Arrow'
#mergeBlocks1Display.GlyphTableIndexArray = 'None'
#mergeBlocks1Display.GaussianRadius = 0.001850000000000001
#mergeBlocks1Display.SetScaleArray = [None, '']
#mergeBlocks1Display.ScaleTransferFunction = 'PiecewiseFunction'
#mergeBlocks1Display.OpacityArray = [None, '']
#mergeBlocks1Display.OpacityTransferFunction = 'PiecewiseFunction'
#mergeBlocks1Display.DataAxesGrid = 'GridAxesRepresentation'
#mergeBlocks1Display.PolarAxes = 'PolarAxesRepresentation'
#mergeBlocks1Display.ScalarOpacityUnitDistance = 0.0065557524973083444

# hide data in view
#Hide(tRACEcgns, renderView1)

# update the view to ensure updated data information
#renderView1.Update()

# save data
SaveData('solution.vtk', proxy=mergeBlocks1, FileType='Ascii')

#### saving camera placements for all active views

# current camera placement for renderView1
#renderView1.CameraPosition = [0.17330032824415897, -0.23223949263701013, 0.25988155242216526]
#renderView1.CameraFocalPoint = [-0.053089408331868204, 0.41895306780119, -0.4032137412963748]
#renderView1.CameraViewUp = [0.047228799639132624, 0.7206988812780291, 0.6916375951386998]
#renderView1.CameraParallelScale = 0.24757491740843654

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
