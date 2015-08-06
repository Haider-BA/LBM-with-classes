#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_
#include <vector>
#include "ParticleNode.hpp"

class Particle {
 public:
  Particle(std::size_t num_nodes
    , double radius
    , double stiffness
    , double center_x
    , double center_y);

  virtual ~Particle() = default;

  void AddNode(double x
    , double y
    , double x_ref
    , double y_ref
    , double u_x
    , double u_y
    , double force_x
    , double force_y);

  ParticleNode center;
  std::vector<ParticleNode> nodes;
};
#endif  // PARTICLE_HPP_
