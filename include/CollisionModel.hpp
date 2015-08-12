#ifndef COLLISION_MODEL_HPP_
#define COLLISION_MODEL_HPP_
#include <vector>
#include "LatticeModel.hpp"

class CollisionModel {
 public:
  /**
   * Constructor: Creates collision model with the same density at each node
   * \param lm lattice model used for simulation
   * \param initial_density initial density of the lattice
   */
  CollisionModel(LatticeModel &lm
    , double initial_density);

  /**
   * Constructor: Creates collision model with the same density at each node
   * \param lm lattice model used for simulation
   * \param initial_density initial density of the lattice
   */
  CollisionModel(LatticeModel &lm
    , const std::vector<double> &initial_density);

  // https://stackoverflow.com/questions/353817/should-every-class-have-a-
  // virtual-destructor
  /**
   * Virtual destructor since we are deriving from this class
   */
  virtual ~CollisionModel()= default;

  /**
   * Calculates equilibrium distribution function according to LBIntro
   */
  void ComputeEq();

  /**
   * Compute density at each node by summing up its distribution functions
   * \param lattice 2D vector containing distribution functions
   * \return density of lattice stored row-wise in a 1D vector
   */
  std::vector<double> ComputeRho(const std::vector<std::vector<double>> &df);

  /**
   * Pure virtual function to compute the macroscopic properties of the lattice
   * depending on the equation, density and velocity for Navier-Stokes, only
   * density for Convection-diffusion equation
   * This is used to unify function calling in LatticeBoltzmann TakeStep()
   * method
   * \param df Particle distribution functions of the lattice stored row-wise
   *        in a 2D vector
   */
  virtual void ComputeMacroscopicProperties(
      const std::vector<std::vector<double>> &df) = 0;

  /**
   * Pure virtual function to compute collision step and apply force step
   * according to "A new scheme for source term in LBGK model for
   * convection–diffusion equation" and Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  virtual void Collide(std::vector<std::vector<double>> &df) = 0;

  /**
   * Adds a node to exclude it from the collision step
   * \param n index of the node in the lattice
   */
  void AddNodeToSkip(std::size_t n);

  /**
   * Equilibrium distribution function stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> edf;

  /**
   * Density stored row-wise in a 1D vector.
   */
  std::vector<double> rho;

  /**
   * Skips the collision step for the node if it is a full-way bounceback node
   */
  std::vector<bool> skip;

 protected:
  /**
   * Lattice model to handle number of rows, columns, dimensions, directions,
   * velocity
   */
  LatticeModel &lm_;

  /**
   * Relaxation time
   */
  double tau_;

  /**
   * Speed of sound in lattice
   */
  double c_;

  /**
   * Square of speed of sound in lattice, used to simplify computations in the
   * collision step
   */
  double cs_sqr_ = c_ * c_ / 3.0;
};
#endif // COLLISION_MODEL_HPP_
