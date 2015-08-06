#include "ParticleNode.hpp"

ParticleNode::ParticleNode(double x
  , double y
  , double x_ref
  , double y_ref)
  : coord {{x, y}},
    coord_ref {{x_ref, y_ref}},
    u {{0.0, 0.0}},
    force {{0.0, 0.0}}
{}

ParticleNode::ParticleNode(double x
  , double y
  , double x_ref
  , double y_ref
  , double u_x
  , double u_y
  , double force_x
  , double force_y)
  : coord {{x, y}},
    coord_ref {{x_ref, y_ref}},
    u {{u_x, u_y}},
    force {{force_x, force_y}}
{}
