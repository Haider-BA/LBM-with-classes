#ifndef IMMERSED_BOUNDARY_METHOD_
#define IMMERSED_BOUNDARY_METHOD_
#include <vector>

#include "Particle.hpp"

class ImmersedBoundaryMethod {
 public:
  ImmersedBoundaryMethod(int interpolation_stencil
    , std::vector<std::vector<double>> &lattice_force
    , std::vector<std::vector<double>> &lattice_velocity);

  ~ImmersedBoundaryMethod() = default;

  void AddParticle(Particle* particle);

  void InterpolateFluidVelocity();

  void SpreadForce();

  std::vector<Particle*> particles;

  std::vector<std::vector<double>> &fluid_force;

  std::vector<std::vector<double>> &fluid_velocity;

 private:
  int interpolation_stencil_;

};
#endif  // IMMERSED_BOUNDARY_METHOD_
