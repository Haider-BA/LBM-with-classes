#ifndef PARTICLE_RIGID_HPP_
#define PARTICLE_RIGID_HPP_
#include "Particle.hpp"

class ParticleRigid: public Particle {
 public:
  ParticleRigid(std::size_t num_nodes
    , double radius
    , double stiffness
    , double center_x
    , double center_y);

  ~ParticleRigid() = default;

  void ComputeForces();

 protected:
};

#endif  // PARTICLE_RIGID_HPP_
