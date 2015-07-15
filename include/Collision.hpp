#ifndef COLLISION_HPP_
#define COLLISION_HPP_
#include <vector>
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"

class Collision {
 public:
  // "Overloading" pure virtual function doesn't work
  // https://stackoverflow.com/questions/15827632/overload-of-pure-virtual-
  // function
  // virtual void InitSource() = 0;
  Collision();
  /**
   * Constructor: Creates collision model with the same density at each node
   * \param lm lattice model used for simulation
   * \param initial_density initial density of the lattice
   */
  Collision(LatticeModel &lm
    , double initial_density);

  /**
   * Constructor: Creates collision model with variable density at each node
   * \param lm lattice model used for simulation
   * \param initial_density initial density of the lattice
   */
  Collision(LatticeModel &lm
    , const std::vector<double> &initial_density);

  // Since we are deriving from this class, need virtual destructor
  // https://stackoverflow.com/questions/353817/should-every-class-have-a-
  // virtual-destructor
  /**
   * Virtual destructor since we are deriving from this class
   */
  virtual ~Collision() = default;

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
   * Pure virtual function for removing source term
   */
  virtual void KillSource() = 0;

  /**
   * Equilibrium distribution function stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> lattice_eq;

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
  bool is_implemented;

 protected:
  // have to use const reference as default constructor will create temporary
  // LatticeModel object
  // https://stackoverflow.com/questions/17905101/invalid-initialization-of-non-
  // const-reference-of-type-stdvectordouble-fro
  LatticeModel &lm_;
  double tau_;
  double c_;
  double cs_sqr_ = c_ * c_ / 3.0;
};
#endif  // COLLISION_HPP_
