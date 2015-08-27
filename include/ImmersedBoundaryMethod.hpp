#ifndef IMMERSED_BOUNDARY_METHOD_
#define IMMERSED_BOUNDARY_METHOD_
#include <vector>
#include "LatticeModel.hpp"
#include "Particle.hpp"

class ImmersedBoundaryMethod {
 public:
  /**
   * Constructor: Creates the immersed boundary method solver with the
   * selected interpolation stencil, reference to lattice force and velocity
   * \param interpolation_stencil the interpolation stencil to be used for
   *        SpreadForce() and InterpolateFluidVelocity() methods
   * \param lattice_force reference to lattice force
   * \param lm reference to lattice model, provides information on number of
   *        rows, columns, space step, time step and velocity of the lattice
   */
  ImmersedBoundaryMethod(int interpolation_stencil
    , std::vector<std::vector<double>> &lattice_force
    , LatticeModel &lm);

  /**
   * Destructor
   */
  ~ImmersedBoundaryMethod() = default;

  /**
   * Adds a particle to the immersed boundary method solver
   * \param particle pointer to the particle to be added
   */
  void AddParticle(Particle *particle);

  /**
   * Spread particle forces to the lattice using the specified interpolation
   * stencil
   * "Introduction to immersed boundary methods"
   */
  void SpreadForce();

  /**
   * Interpolate the fluid velocity to the boundary nodes in the particles using
   * the specified interpolation stencil
   * "Introduction to immersed boundary methods"
   */
  void InterpolateFluidVelocity();

  /**
   * Update node positions of the particles based on node velocities
   * "Introduction to immersed boundary methods"
   */
  void UpdateParticlePosition();

  std::vector<Particle*> particles;

  std::vector<std::vector<double>> &fluid_force;

 private:
  int interpolation_stencil_;

  LatticeModel &lm_;

};
#endif  // IMMERSED_BOUNDARY_METHOD_
