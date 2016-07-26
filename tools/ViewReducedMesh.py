try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

ReducedMesh_xmf = GetActiveSource()
ExtractSurface1 = ExtractSurface()

RenderView1 = GetRenderView()
DataRepresentation1 = GetDisplayProperties(ReducedMesh_xmf)
DataRepresentation2 = Show()
DataRepresentation2.ColorArrayName = 'reduced mesh'
DataRepresentation2.ScaleFactor = 0.24
DataRepresentation2.SelectionPointFieldDataArrayName = 'reduced mesh'
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Threshold1 = Threshold()

a1_reducedmesh_PVLookupTable = GetLookupTableForArray( "reduced mesh", 1 )

RenderView1.CameraViewUp = [-0.24535093172717765, 0.9533130041025328, -0.1760603206563314]

DataRepresentation1.Visibility = 0

DataRepresentation2.LookupTable = a1_reducedmesh_PVLookupTable

Threshold1.Scalars = ['POINTS', 'reduced mesh']
Threshold1.ThresholdRange = [0.0, 2.0]

Threshold1.ThresholdRange = [2.0, 2.0]

DataRepresentation3 = Show()
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation3.SelectionPointFieldDataArrayName = 'reduced mesh'
DataRepresentation3.ColorArrayName = 'reduced mesh'
DataRepresentation3.ScalarOpacityUnitDistance = 0.33727814214971275
DataRepresentation3.ScaleFactor = 0.04958867281675339

SetActiveSource(ReducedMesh_xmf)
Threshold2 = Threshold()

RenderView1.CameraClippingRange = [7.128036687197317, 8.180663808067102]

DataRepresentation2.Visibility = 0

DataRepresentation3.Representation = 'Surface With Edges'
DataRepresentation3.ScalarOpacityFunction = []
DataRepresentation3.LookupTable = a1_reducedmesh_PVLookupTable

Threshold2.Scalars = ['POINTS', 'reduced mesh']
Threshold2.ThresholdRange = [0.0, 2.0]

Threshold2.ThresholdRange = [1.0, 1.0]

DataRepresentation4 = Show()
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation4.SelectionPointFieldDataArrayName = 'reduced mesh'
DataRepresentation4.ColorArrayName = 'reduced mesh'
DataRepresentation4.ScalarOpacityUnitDistance = 0.8174976848755081
DataRepresentation4.ScaleFactor = 0.24000000953674316

RenderView1.CameraViewUp = [-0.2975564271576019, 0.941405886694574, 0.15879272386984836]
RenderView1.CameraPosition = [6.6003450468770755, 2.7080706335376585, -3.686683694151925]
RenderView1.CameraClippingRange = [4.085973509822269, 13.01534395974289]

DataRepresentation1.Opacity = 0.4
DataRepresentation1.ColorArrayName = ''
DataRepresentation1.Visibility = 1

DataRepresentation4.Representation = 'Surface With Edges'
DataRepresentation4.ScalarOpacityFunction = []
DataRepresentation4.LookupTable = a1_reducedmesh_PVLookupTable

Render()
