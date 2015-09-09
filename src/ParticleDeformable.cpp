#include "ParticleDeformable.hpp"
#include "LatticeModel.hpp"

ParticleDeformable::ParticleDeformable(double stiffness
  , double bending
  , std::size_t num_nodes
  , double center_x
  , double center_y
  , LatticeModel &lm)
  : Particle(stiffness, num_nodes, center_x, center_y, lm)
{}

void ParticleDeformable::ComputeForces()
{
  ;
}
