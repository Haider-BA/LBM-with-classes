#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#### function to set up stream tracers and rename them
def SetupTrace(obj, num_nodes_x, res, x_pos, y_pos1, y_pos2, name):
    obj.Vectors = ['POINTS', 'velocity_vector']
    obj.InterpolatorType = 'Interpolator with Cell Locator'
    obj.InitialStepLength = 0.01
    obj.MaximumStepLength = 0.01
    obj.MaximumSteps = 80000
    streamTracer1.MaximumStreamlineLength = num_nodes_x
    obj.SeedType.Point1 = [x_pos, y_pos1, 0.0]
    obj.SeedType.Point2 = [x_pos, y_pos2, 0.0]
    obj.SeedType.Resolution = res
    RenameSource(name, obj)
    pass

#### simulation parameters
nx = 800  # length of channel
ny = 200  # height of channel
radius = ny / 5  # cylinder radius
l_edge = 0  # left edge position
w_cylinder = nx / 4 - radius * 1.25  # west of cylinder
e_cylinder = nx / 4 + radius * 1.25  # east of cylinder
n_cylinder = ny / 2 + radius  # north of cylinder (exact)
s_cylinder = ny / 2 - radius  # south of cylinder (exact)
d_cylinder = nx / 4 + radius * 2  # downstream of cylinder

# find source
fluid_t = FindSource('fluid_t*')
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# get display properties
fluid_tDisplay = GetDisplayProperties(fluid_t, view=renderView1)
# set scalar coloring
ColorBy(fluid_tDisplay, ('POINTS', 'velocity_vector'))
# get color transfer function/color map for 'velocityvector'
velocityvectorLUT = GetColorTransferFunction('velocityvector')
# get opacity transfer function/opacity map for 'velocityvector'
velocityvectorPWF = GetOpacityTransferFunction('velocityvector')
# Apply a preset using its name. Note this may not work as expected when presets
# have duplicate names.
velocityvectorLUT.ApplyPreset('Blue to Red Rainbow', True)
# rescale color and/or opacity maps used to exactly fit the current data range
fluid_tDisplay.RescaleTransferFunctionToDataRange(False)

# set active source
SetActiveSource(fluid_t)
# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=fluid_t,
    SeedType='High Resolution Line Source')
SetupTrace(streamTracer1,
    1.5 * nx,
    20,
    l_edge,
    0,
    ny - 1,
    'Left edge')
# show data in view
streamTracer1Display = Show(streamTracer1, renderView1)
# trace defaults for the display properties.
streamTracer1Display.ColorArrayName = [None, '']

# set active source
SetActiveSource(fluid_t)
# create a new 'Stream Tracer'
streamTracer2 = StreamTracer(Input=fluid_t,
    SeedType='High Resolution Line Source')
SetupTrace(streamTracer2,
    1.5 * nx,
    10,
    w_cylinder,
    s_cylinder,
    n_cylinder,
    'West of cylinder')
# show data in view
streamTracer2Display = Show(streamTracer2, renderView1)
# trace defaults for the display properties.
streamTracer2Display.ColorArrayName = [None, '']

# set active source
SetActiveSource(fluid_t)
# create a new 'Stream Tracer'
streamTracer3 = StreamTracer(Input=fluid_t,
    SeedType='High Resolution Line Source')
SetupTrace(streamTracer3,
    1.5 * nx,
    10,
    e_cylinder,
    s_cylinder,
    n_cylinder,
    'East of cylinder')
# show data in view
streamTracer3Display = Show(streamTracer3, renderView1)
# trace defaults for the display properties.
streamTracer3Display.ColorArrayName = [None, '']

# set active source
SetActiveSource(fluid_t)
# create a new 'Stream Tracer'
streamTracer4 = StreamTracer(Input=fluid_t,
    SeedType='High Resolution Line Source')
SetupTrace(streamTracer4,
    1.5 * nx,
    10,
    d_cylinder,
    s_cylinder,
    n_cylinder,
    'Downstream of cylinder')
# show data in view
streamTracer4Display = Show(streamTracer4, renderView1)
# trace defaults for the display properties.
streamTracer4Display.ColorArrayName = [None, '']

# create a new 'Cylinder'
cylinder1 = Cylinder()
# Properties modified on cylinder1
cylinder1.Resolution = 360
cylinder1.Radius = radius
cylinder1.Center = [nx / 4 - 1, ny / 2 - 1, 0.0]
# create a new 'Transform'
transform1 = Transform(Input=cylinder1)
transform1.Transform = 'Transform'
# Properties modified on transform1.Transform
transform1.Transform.Translate = [0.0, ny / 2 - 1, 0.0]
transform1.Transform.Rotate = [90.0, 0.0, 0.0]
# show data in view
transform1Display = Show(transform1, renderView1)
# set rotated cylinder to be black
transform1Display.DiffuseColor = [0.0, 0.0, 0.0]

#### saving camera placements for all active views
# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [398.2520875930786, 99.49999713897705,
        1597.7643741359389]
renderView1.CameraFocalPoint = [398.2520875930786, 99.49999713897705, 0.0]
renderView1.CameraParallelScale = 413.5318496126904
