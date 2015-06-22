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

  /**
   * Returns density lattice
   * \return density lattice
   */
  std::vector<double> GetRho() const;

  /**
   * Returns velocity lattice
   * \return velocity lattice
   */
  std::vector<std::vector<double>> GetVelocity() const;

  /** \brief
   *
   * \return void
   *
   */
  void ComputeEq();

  /** \brief
   *
   * \return void
   *
   */
  void Collide(std::vector<std::vector<double>> &lattice);

  double InnerProduct(const std::vector<double> &a
    , const std::vector<double> &b);

  /** \brief
   *
   * \return virtual void
   *
   */
  virtual void ApplyForce(std::vector<std::vector<double>> &lattice) = 0;

  /**
   * Equilibrium distribution function stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> lattice_eq;
 protected:
  LatticeModel &lm_;

  /**
   * Lattice velocity stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> u_;

  /**
   * Density stored row-wise in a 1D vector.
   */
  std::vector<double> rho_;
  double tau_;
  double c_;
  double cs_sqr_;
};
#endif  // COLLISION_HPP_
