#include "ParticleRigid.hpp"

ParticleRigid::ParticleRigid(double stiffness
  , std::size_t num_nodes
  , double center_x
  , double center_y)
  : Particle(stiffness, num_nodes, center_x, center_y)
{}

void ParticleRigid::ComputeForces()
{}
