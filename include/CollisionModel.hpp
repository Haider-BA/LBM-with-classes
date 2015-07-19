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
   * Pure virtual function to compute collision step and apply force step
   * according to "A new scheme for source term in LBGK model for
   * convection–diffusion equation" and Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  virtual void Collide(std::vector<std::vector<double>> &lattice) = 0;

  /**
   * Equilibrium distribution function stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> edf;

  /**
   * Density stored row-wise in a 1D vector.
   */
  std::vector<double> rho;

 protected:
  LatticeModel &lm_;
  double tau_;
  double c_;
  double cs_sqr_ = c_ * c_ / 3.0;
};
#endif // COLLISION_MODEL_HPP_
