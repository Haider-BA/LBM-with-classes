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
   * Immersed Boundary Method Interpolation Stencil
   * http://lbmworkshop.com/wp-content/uploads/2011/09/2011-08-25_Edmonton_IBM.pdf
   * \param x position in one of the lattice directions
   * \return contribution towards the interpolation function in one of the lattice
   *         directions
   */
  double Phi2(double x);

  /**
   * Immersed Boundary Method Interpolation Stencil
   * http://lbmworkshop.com/wp-content/uploads/2011/09/2011-08-25_Edmonton_IBM.pdf
   * \param x position in one of the lattice directions
   * \return contribution towards the interpolation function in one of the lattice
   *         directions
   */
  double Phi3(double x);

  /**
   * Immersed Boundary Method Interpolation Stencil
   * http://lbmworkshop.com/wp-content/uploads/2011/09/2011-08-25_Edmonton_IBM.pdf
   * \param x position in one of the lattice directions
   * \return contribution towards the interpolation function in one of the lattice
   *         directions
   */
  double Phi4(double x);

  /**
   * Approximates Dirac delta-function of IBM using Phi2
   * \param x x-position in lattice
   * \param y y-position in lattice
   * \return approximated Dirac delta-function value
   */
  double Dirac2(double x, double y);

  /**
   * Approximates Dirac delta-function of IBM using Phi3
   * \param x x-position in lattice
   * \param y y-position in lattice
   * \return approximated Dirac delta-function value
   */
  double Dirac3(double x, double y);

  /**
   * Approximates Dirac delta-function of IBM using Phi4
   * \param x x-position in lattice
   * \param y y-position in lattice
   * \return approximated Dirac delta-function value
   */
  double Dirac4(double x, double y);

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
