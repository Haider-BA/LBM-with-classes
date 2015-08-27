#include "ParticleRigid.hpp"
#include <iostream>

ParticleRigid::ParticleRigid(double stiffness
  , std::size_t num_nodes
  , double center_x
  , double center_y)
  : Particle(stiffness, num_nodes, center_x, center_y)
{}

void ParticleRigid::ComputeForces()
{
  auto area = area_ / nodes.size();
  for (auto &node : nodes) {
    node.force = {0.0, 0.0};
    node.force[0] = -stiffness_ * area * (node.coord[0] - node.coord_ref[0]);
    node.force[1] = -stiffness_ * area * (node.coord[1] - node.coord_ref[1]);
  }  // node
}
