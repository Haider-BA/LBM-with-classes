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
   * Compute density at each node by summing up its distribution functions
   * \param lattice 2D vector containing distribution functions
   * \return density of lattice stored row-wise in a 1D vector
   */
  std::vector<double> ComputeRho(const std::vector<std::vector<double>> &df);

  virtual void ComputeMacroscopicProperties(
      const std::vector<std::vector<double>> &df) = 0;

  /**
   * Pure virtual function to compute collision step and apply force step
   * according to "A new scheme for source term in LBGK model for
   * convection–diffusion equation" and Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  virtual void Collide(std::vector<std::vector<double>> &lattice) = 0;

  /** \brief
   *
   * \param n std::size_t
   * \return void
   *
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

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  std::vector<bool> skip;

 protected:
  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  LatticeModel &lm_;

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  double tau_;

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  double c_;

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  double cs_sqr_ = c_ * c_ / 3.0;
};
#endif // COLLISION_MODEL_HPP_
