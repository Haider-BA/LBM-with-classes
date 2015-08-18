#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_
#include <vector>
#include "ParticleNode.hpp"

class Particle {
 public:
  Particle(double stiffness
    , std::size_t num_nodes
    , double center_x
    , double center_y);

  virtual ~Particle() = default;

  std::size_t GetNumberOfNodes() const;

  void AddNode(double x
    , double y
    , double x_ref
    , double y_ref
    , double u_x
    , double u_y
    , double force_x
    , double force_y);

  virtual void ComputeForces() = 0;

  void CreateCylinder(double radius);

  ParticleNode center;

  std::vector<ParticleNode> nodes;
 protected:
  double pi_ = 3.14159265;
  std::size_t number_of_nodes_;
};
#endif  // PARTICLE_HPP_
