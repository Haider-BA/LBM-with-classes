#include "ParticleRigid.hpp"
#include <iostream>
#include "LatticeModel.hpp"

ParticleRigid::ParticleRigid(double stiffness
  , std::size_t num_nodes
  , double center_x
  , double center_y
  , LatticeModel &lm)
  : Particle(stiffness, num_nodes, center_x, center_y, lm)
{}

void ParticleRigid::ComputeForces()
{
  // change to store coordinates on particle as real coordinates
  const auto area = area_ / nodes.size();
  for (auto &node : nodes) {
    node.force[0] = -stiffness_ * area * (node.coord[0] - node.coord_ref[0]);
    node.force[1] = -stiffness_ * area * (node.coord[1] - node.coord_ref[1]);
  }  // node
}
