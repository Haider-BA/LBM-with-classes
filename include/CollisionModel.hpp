#ifndef COLLISIONMODEL_HPP_
#define COLLISIONMODEL_HPP_
#include <vector>
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"

class CollisionModel {
 public:
  // "Overloading" pure virtual function doesn't work
  // https://stackoverflow.com/questions/15827632/overload-of-pure-virtual-
  // function
  // virtual void InitSource() = 0;

  /**
   * Constructor: Creates collision model with the same density at each node
   * \param lm lattice model used for simulation
   * \param initial_density initial density of the lattice
   */
  CollisionModel(LatticeModel &lm);

  /**
   * Constructor: Creates collision model with variable density at each node
   * \param lm lattice model used for simulation
   * \param initial_density initial density of the lattice
   */
  CollisionModel(LatticeModel &lm
    , const std::vector<double> &initial_density);

  // Since we are deriving from this class, need virtual destructor
  // https://stackoverflow.com/questions/353817/should-every-class-have-a-
  // virtual-destructor
  /**
   * Virtual destructor since we are deriving from this class
   */
  virtual ~CollisionModel() = default;

  /**
   * Calculates equilibrium distribution function according to LBIntro
   */
  void ComputeEq(std::vector<std::vector<double>> &lattice_eq
    , const std::vector<double> &rho);

  /**
   * Virtual function to compute collision step and apply force step depending
   * on the model used according to Guo2002
   * \param lattice 2D vector containing distribution functions
   */
//  virtual void CollideNS(std::vector<std::vector<double>> &lattice);

  /**
   * Virtual function to compute collision step and apply force step
   * according to "A new scheme for source term in LBGK model for
   * convection–diffusion equation"
   * \param lattice 2D vector containing distribution functions
   */
//  virtual void CollideCD(std::vector<std::vector<double>> &lattice) {};

  /**
   * Equilibrium distribution function stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> f_eq;

  std::vector<std::vector<double>> g_eq;

  /**
   * Density stored row-wise in a 1D vector.
   */
  std::vector<double> rho_f;

  std::vector<double> rho_g;

  /**
   * Indicates which collision model is implemented
   */
  bool is_ns;
  bool is_cd;

 protected:
  // have to use const reference as default constructor will create temporary
  // LatticeModel object
  // https://stackoverflow.com/questions/17905101/invalid-initialization-of-non-
  // const-reference-of-type-stdvectordouble-fro
  // no longer using const reference as using temporary will give warning
  /**
   * Lattice model which contains the number of dimensions, directions, rows,
   * columns, and the discrete velocity vectors
   */
  LatticeModel &lm_;

  /**
   * Relaxation-time
   */
//  double tau_;

  /**
   * Lattice speed (dx / dt)
   */
  double c_;

  /**
   * Square of speed of sound in lattice
   */
  double cs_sqr_ = c_ * c_ / 3.0;
};
#endif  // COLLISIONMODEL_HPP_
