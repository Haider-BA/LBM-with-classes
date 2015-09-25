$first = 0;
$last = 500;
$step = 1;

sub anim{
  gfx read elem lbm$first;
  for( $a = $step; $a <= $last; $a = $a + $step ){
    gfx read elem "lbm$a";
    gfx update;
  }
}

sub animprint{
  gfx read elem lbm$first;
  gfx print file lbm0.png anti 8 force width 720 height 480;
  for( $a = $step; $a <= $last; $a = $a + $step ){
    gfx read elem "lbm$a";
    $b = $a/$step;
    gfx update;
    gfx print file "lbm$b.png" anti 8 force  width 720 height 480;
  }
}

gfx read node lbm;
gfx read elem lbm$first;

gfx define field velocity component u_x u_y;
gfx define field vmag magnitude field velocity;

gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout 2d ortho_axes z -y eye_spacing 0.25 width 720 height 480;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point 14.5 9.5 70.723 interest_point 14.5 9.5 0 up_vector 0 -1 -0 view_angle 24.7183 near_clipping_plane 0.70723 far_clipping_plane 252.74 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0 0.1 extend_above extend_below rainbow colour_range 0 1 component 1;
# gfx create spectrum "vectors" clear overwrite_colour;
# gfx modify spectrum "vectors" linear range 0 1 extend_above extend_below monochrome colour_range 0 1 component 1;
gfx create spectrum "streamlines" clear overwrite_colour;
gfx modify spectrum "streamlines" linear range 0 0 extend_above extend_below alpha colour_range 0 1 component 1;


# gfx create spectrum "porosity";
# gfx modify spectrum "porosity" clear overlay_colour;
# gfx modify spectrum "porosity" linear range 0 1 extend_above extend_below monochrome colour_range 0 1 component 1;
# gfx modify spectrum "porosity" log exaggeration 100 right range 0 1 extend_above extend_below alpha colour_range 0 0.3 component 1;

gfx modify g_element "/" general clear;
gfx modify g_element "/" lines as "border" coordinate coordinates tessellation default LOCAL select_on material default selected_material default_selected;

$x_pos = 0;
$num_lines = 20;
# slightly smaller spacing so all lines can be displayed
$spacing = 1 / ($num_lines + 1);
for ($i = 0; $i < $num_lines; $i = $i + 1) {
  $y_pos = ($i + 1) * $spacing;
  gfx modify g_element "/" streamlines as "vel line $i" coordinate coordinates discretization "1*1*1" tessellation NONE LOCAL exact_xi native_discretization pressure xi $x_pos,$y_pos,0.5 line vector velocity length 100 width 2 travel_scalar select_on material black selected_material black spectrum "streamlines";
}

# gfx modify g_element "/" streamlines as "velocity vectors" coordinate coordinates discretization "1*1*1" tessellation NONE LOCAL cell_corners native_discretization pressure line vector velocity length 10 width 2 travel_scalar select_on material black selected_material default_selected spectrum "vectors";
gfx modify g_element "/" surfaces as "velocity magnitude" coordinate coordinates tessellation default LOCAL native_discretization pressure select_on material default data vmag spectrum default selected_material default_selected render_shaded;

# gfx modify g_element "/" element_points coordinate coordinates discretization "1*1*1" tessellation NONE LOCAL glyph cube_solid general size "1*1*1" centre 0,0,0 font default use_elements cell_centres native_discretization porosity select_on material default data porosity spectrum "porosity" selected_material default_selected;

gfx modify spectrum default autorange;

# gfx modify g_element "/" general clear;
# gfx modify g_element "/" lines coordinate coordinates tessellation default LOCAL select_on material default selected_material default_selected;
# gfx modify g_element "/" surfaces coordinate coordinates tessellation default LOCAL native_discretization pressure select_on material default data u_x spectrum default selected_material default_selected render_shaded;

# gfx modify g_element "/" general clear;
# gfx modify g_element "/" lines coordinate coordinates tessellation default LOCAL select_on material default selected_material default_selected;
# gfx modify g_element "/" streamlines coordinate coordinates discretization "30*20*20" tessellation NONE LOCAL cell_random line vector velocity length 10 width 2 no_data select_on material black selected_material default_selected;
# gfx modify g_element "/" surfaces coordinate coordinates tessellation default LOCAL native_discretization pressure select_on material default data vmag spectrum default selected_material default_selected render_shaded;

# avconv -framerate 10 -f image2 -i lbm%d.png -c:v h264 -crf 1 out.mov
