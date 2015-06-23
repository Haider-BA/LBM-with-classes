#ifndef COLLISION_HPP_
#define COLLISION_HPP_
#include <vector>
#include "LatticeModel.hpp"

class Collision {
 public:
  // "Overloading" pure virtual function doesn't work
  // https://stackoverflow.com/questions/15827632/overload-of-pure-virtual-
  // function
  //virtual void InitSource() = 0;
  Collision(LatticeModel &lm
    , double initial_density
    , const std::vector<double> &initial_velocity);

  ~Collision() = default;

  /**
   * Calculates equilibrium distribution function according to LBIntro
   */
  void ComputeEq();

  /**
   * Computes collision step according to LBIntro
   * \param lattice 2D vector containing distribution functions
   */
  void Collide(std::vector<std::vector<double>> &lattice);

  /**
   * Pure virtual function for forcing term
   * \param lattice 2D vector containing distribution functions
   */
  virtual void ApplyForce(std::vector<std::vector<double>> &lattice) = 0;

  /**
   * Does dot product between 2 vectors of equal length
   * \param a first vector
   * \param b second vector
   * \return dot product of a and b
   */
  double InnerProduct(const std::vector<double> &a
  , const std::vector<double> &b);

  /**
   * Equilibrium distribution function stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> lattice_eq;

  /**
   * Lattice velocity stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> u;

  /**
   * Density stored row-wise in a 1D vector.
   */
  std::vector<double> rho;

 protected:
  LatticeModel &lm_;
  double tau_;
  double c_;
  double cs_sqr_;
};
#endif  // COLLISION_HPP_
