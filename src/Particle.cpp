#include "Particle.hpp"
#include <cmath>
#include "ParticleNode.hpp"

Particle::Particle(double stiffness
  , std::size_t num_nodes
  , double center_x
  , double center_y)
  : center {ParticleNode(center_x, center_y, center_x, center_y)},
    nodes {},
    area_ {0.0},
    stiffness_ {stiffness},
    number_of_nodes_ {num_nodes}
{}

std::size_t Particle::GetNumberOfNodes() const
{
  return number_of_nodes_;
}

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

void Particle::CreateCylinder(double radius)
{
  for (auto i = 0u; i < number_of_nodes_; ++i) {
    auto x = center.coord[0] + radius * sin(2.0 * pi_ * static_cast<double>(i) /
        number_of_nodes_);
    auto y = center.coord[1] + radius * cos(2.0 * pi_ * static_cast<double>(i) /
        number_of_nodes_);
    Particle::AddNode(x, y, x, y, 0.0, 0.0, 0.0, 0.0);
    area_ = 2 * pi_ * radius;
  }
}
