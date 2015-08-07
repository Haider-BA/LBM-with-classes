#include "ParticleRigid.hpp"

ParticleRigid::ParticleRigid(std::size_t num_nodes
  , double radius
  , double stiffness
  , double center_x
  , double center_y)
  : Particle(num_nodes, radius, stiffness, center_x, center_y)
{}

void ParticleRigid::ComputeForces()
{}
