#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
fluid_t = FindSource('fluid_t*')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# get color transfer function/color map for 'densitydifference'
densitydifferenceLUT = GetColorTransferFunction('densitydifference')

# get opacity transfer function/opacity map for 'densitydifference'
densitydifferencePWF = GetOpacityTransferFunction('densitydifference')

# set active source
SetActiveSource(fluid_t)
# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=fluid_t,
    SeedType='High Resolution Line Source')
streamTracer1.Vectors = ['POINTS', 'velocity_vector']
streamTracer1.MaximumStreamlineLength = 799.0
streamTracer1.SeedType.Point1 = [0.0, 0.0, 0.0]
streamTracer1.SeedType.Point2 = [0.0, 199.0, 0.0]
streamTracer1.SeedType.Resolution = 20
streamTracer1.InterpolatorType = 'Interpolator with Cell Locator'
streamTracer1.InitialStepLength = 0.01
streamTracer1.MaximumStepLength = 0.01
streamTracer1.MaximumSteps = 80000
RenameSource('Left edge', streamTracer1)
# show data in view
streamTracer1Display = Show(streamTracer1, renderView1)
# trace defaults for the display properties.
streamTracer1Display.ColorArrayName = [None, '']

# set active source
SetActiveSource(fluid_t)
# create a new 'Stream Tracer'
streamTracer2 = StreamTracer(Input=fluid_t,
    SeedType='High Resolution Line Source')
streamTracer2.Vectors = ['POINTS', 'velocity_vector']
streamTracer2.MaximumStreamlineLength = 799.0
streamTracer2.SeedType.Point1 = [150.0, 59.0, 0.0]
streamTracer2.SeedType.Point2 = [150.0, 140.0, 0.0]
streamTracer2.SeedType.Resolution = 10
streamTracer2.InterpolatorType = 'Interpolator with Cell Locator'
streamTracer2.InitialStepLength = 0.01
streamTracer2.MaximumStepLength = 0.01
streamTracer2.MaximumSteps = 80000
RenameSource('Left cylinder', streamTracer2)
# show data in view
streamTracer2Display = Show(streamTracer2, renderView1)
# trace defaults for the display properties.
streamTracer2Display.ColorArrayName = [None, '']

# set active source
SetActiveSource(fluid_t)
# create a new 'Stream Tracer'
streamTracer3 = StreamTracer(Input=fluid_t,
    SeedType='High Resolution Line Source')
streamTracer3.Vectors = ['POINTS', 'velocity_vector']
streamTracer3.MaximumStreamlineLength = 799.0
streamTracer3.SeedType.Point1 = [250.0, 59.0, 0.0]
streamTracer3.SeedType.Point2 = [250.0, 140.0, 0.0]
streamTracer3.SeedType.Resolution = 10
streamTracer3.InterpolatorType = 'Interpolator with Cell Locator'
streamTracer3.InitialStepLength = 0.01
streamTracer3.MaximumStepLength = 0.01
streamTracer3.MaximumSteps = 80000
RenameSource('Right cylinder', streamTracer3)
# show data in view
streamTracer3Display = Show(streamTracer3, renderView1)
# trace defaults for the display properties.
streamTracer3Display.ColorArrayName = [None, '']

# set active source
SetActiveSource(fluid_t)
# create a new 'Stream Tracer'
streamTracer4 = StreamTracer(Input=fluid_t,
    SeedType='High Resolution Line Source')
streamTracer4.Vectors = ['POINTS', 'velocity_vector']
streamTracer4.MaximumStreamlineLength = 799.0
streamTracer4.SeedType.Point1 = [350.0, 59.0, 0.0]
streamTracer4.SeedType.Point2 = [350.0, 140.0, 0.0]
streamTracer4.SeedType.Resolution = 10
streamTracer4.InterpolatorType = 'Interpolator with Cell Locator'
streamTracer4.InitialStepLength = 0.01
streamTracer4.MaximumStepLength = 0.01
streamTracer4.MaximumSteps = 80000
RenameSource('Downstream cylinder', streamTracer4)
# show data in view
streamTracer4Display = Show(streamTracer4, renderView1)
# trace defaults for the display properties.
streamTracer4Display.ColorArrayName = [None, '']

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [398.2520875930786, 99.49999713897705, 1597.7643741359389]
renderView1.CameraFocalPoint = [398.2520875930786, 99.49999713897705, 0.0]
renderView1.CameraParallelScale = 413.5318496126904

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
