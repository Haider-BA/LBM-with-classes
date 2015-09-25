#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
fluid_t = FindSource('fluid_t*')

# create a new 'Glyph'
glyph1 = Glyph(Input=fluid_t,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'density_difference']
glyph1.Vectors = ['POINTS', 'velocity_vector']
glyph1.ScaleFactor = 21.900000000000002
glyph1.GlyphTransform = 'Transform2'

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
    SeedType='Point Source')
streamTracer1.Vectors = ['POINTS', 'velocity_vector']
streamTracer1.MaximumStreamlineLength = 219.0
streamTracer1.SeedType.Center = [1.0, 30.0, 0.0]
streamTracer1.SeedType.Radius = 0.0
streamTracer1.InterpolatorType = 'Interpolator with Cell Locator'

# show data in view
streamTracer1Display = Show(streamTracer1, renderView1)
# trace defaults for the display properties.
streamTracer1Display.ColorArrayName = [None, '']

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [110.0, 30.0, 10000.0]
renderView1.CameraFocalPoint = [110.0, 30.0, 0.0]
renderView1.CameraParallelScale = 113.40414454507383

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
