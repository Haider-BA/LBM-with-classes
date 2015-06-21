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
  Collision(LatticeModel &lm);

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
  void Collide();

  /** \brief
   *
   * \return virtual void
   *
   */
  virtual void ApplyForce() = 0;
 protected:
  LatticeModel &lm_;
  std::vector<std::vector<double>> lattice_eq_;
  double tau_;
  double c_;
  double cs_sqr_;
};
#endif  // COLLISION_HPP_
