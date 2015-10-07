#include "Particle.hpp"
#include <cmath>
#include "ParticleNode.hpp"

Particle::Particle(double stiffness
  , std::size_t num_nodes
  , double center_x
  , double center_y
  , bool mobility
  , LatticeModel &lm)
  : center {ParticleNode(center_x, center_y, center_x, center_y)},
    nodes {},
    is_mobile {mobility},
    area_ {0.0},
    char_length_ {0.0},
    stiffness_ {stiffness},
    number_of_nodes_ {num_nodes},
    lm_ (lm)
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

void Particle::CreateCylinder(double radius)
{
  // particle coordinates stored in real units
  for (auto i = 0u; i < number_of_nodes_; ++i) {
    const auto x = center.coord[0] + radius * sin(2.0 * pi_ *
        static_cast<double>(i) / number_of_nodes_);
    const auto y = center.coord[1] + radius * cos(2.0 * pi_ *
        static_cast<double>(i) / number_of_nodes_);
    Particle::AddNode(x, y, x, y, 0.0, 0.0, 0.0, 0.0);
  }  // i
  area_ = 2.0 * pi_ * radius;
  char_length_ = radius;
}

void Particle::UpdateReferencePosition()
{
  center.coord_ref = center.coord;
  for (auto i = 0u; i < number_of_nodes_; ++i) {
    nodes[i].coord_ref[0] = center.coord_ref[0] + char_length_ * sin(2.0 * pi_ *
        static_cast<double>(i) / number_of_nodes_);
    nodes[i].coord_ref[1] = center.coord_ref[1] + char_length_ * cos(2.0 * pi_ *
        static_cast<double>(i) / number_of_nodes_);
  }  // i
}

void Particle::ChangeMobility(bool mobility)
{
  is_mobile = mobility;
}

std::size_t Particle::GetNumberOfNodes() const
{
  return number_of_nodes_;
}
