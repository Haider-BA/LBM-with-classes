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

# To make and play a movie, issue these commands into a terminal
# avconv -framerate 10 -f image2 -i lbm%d.png -c:v h264 -crf 1 out.mov
# mplayer out.mov

gfx read node lbm;
gfx read elem lbm$first;

# gfx define field velocity component u_x u_y;
# gfx define field vmag magnitude field velocity;

gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout 2d ortho_axes z -y eye_spacing 0.25 width 720 height 480;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point 14.5 9.5 70.723 interest_point 14.5 9.5 0 up_vector 0 -1 -0 view_angle 24.7183 near_clipping_plane 0.70723 far_clipping_plane 252.74 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

gfx modify g_element "/" general clear;
gfx modify g_element "/" surfaces coordinate coordinates tessellation default LOCAL native_discretization solute select_on material default data solute spectrum default selected_material default_selected render_shaded;

gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0.95 1.05 extend_above extend_below rainbow colour_range 0 1 component 1;

