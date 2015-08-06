#include "Particle.hpp"
#include "ParticleNode.hpp"

Particle::Particle(std::size_t num_nodes
  , double radius
  , double stiffness
  , double center_x
  , double center_y)
  : center {ParticleNode(center_x, center_y, center_x, center_y)},
    nodes {}
{}

void Particle::AddNode(double x
  , double y
  , double x_ref
  , double y_ref
  , double u_x
  , double u_y
  , double force_x
  , double force_y)
{
  nodes.push_back(ParticleNode(x, y, x_ref, y_ref, u_x, u_y, force_x, force_y));
}
