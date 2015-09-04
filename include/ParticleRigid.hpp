#ifndef PARTICLE_RIGID_HPP_
#define PARTICLE_RIGID_HPP_
#include "LatticeModel.hpp"
#include "Particle.hpp"

class ParticleRigid: public Particle {
 public:
  /**
   * Constructor: Creates an rigid particle. Particle walls are quasi-rigid,
   * i.e., wall deformation is small
   * \param stiffness stiffness of the walls, used to calculate restorative
   *        force
   * \param num_nodes number of boundary nodes in the particle
   * \param center_x x-coordinate for particle center
   * \param center_y y-coordinate for particle center
   * \param lm reference to LatticeModel to provide information on number of
   *        rows, columns, dimensions, discrete directions and lattice velocity
   */
  ParticleRigid(double stiffness
    , std::size_t num_nodes
    , double center_x
    , double center_y
    , LatticeModel &lm);

  /**
   * Destructor
   */
  ~ParticleRigid() = default;

  /**
   * Compute particle forces on fluid.
   * Current implementation assumes particles are stationary
   * The immersed boundary-lattice Boltzmann method for solving fluid–particles
   * interaction problems
   */
  void ComputeForces();

  // need to update node reference coordinates through rigid body motion

 protected:
};

#endif  // PARTICLE_RIGID_HPP_
