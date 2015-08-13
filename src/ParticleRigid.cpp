#include "ParticleRigid.hpp"

ParticleRigid::ParticleRigid(double stiffness
  , double center_x
  , double center_y)
  : Particle(stiffness, center_x, center_y)
{}

void ParticleRigid::ComputeForces()
{}
