#ifndef PARTICLE_DEFORMABLE_HPP_
#define PARTICLE_DEFORMABLE_HPP_
#include "LatticeModel.hpp"
#include "Particle.hpp"

class ParticleDeformable: public Particle {
 public:
  ParticleDeformable(double stiffness
    , double bending
    , std::size_t num_nodes
    , double center_x
    , double center_y
    , LatticeModel &lm);

  /**
   * Destructor
   */
  ~ParticleDeformable() = default;

  void ComputeForces();
};

#endif  // PARTICLE_DEFORMABLE_HPP_
