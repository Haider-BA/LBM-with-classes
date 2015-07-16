#ifndef LATTICEBOLTZMANN_HPP_
#define LATTICEBOLTZMANN_HPP_
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeModel.hpp"

class LatticeBoltzmann {
 public:
  /**
   * Constructor: Creates lattice Boltzmann simulation
   * \param t_total total duration of the simulation
   * \param lm lattice model for the simulation
   * \param ns collision model for Navier-Stokes equation
   * \param cd collision model for Convection-Diffusion equation
   */
  LatticeBoltzmann(double t_total
    , LatticeModel &lm
    , CollisionNS &ns
    , CollisionCD &cd);

  /**
   * Destructor
   */
  ~LatticeBoltzmann() = default;

  /**
   * Copy constructor
   */
  LatticeBoltzmann(const LatticeBoltzmann&) = default;

  /**
   * Streams the distribution functions according to LBIntro.
   * \param lattice pre-stream distribution functions
   * \return post-stream distribution functions
   */
  std::vector<std::vector<double>> Stream(
    const std::vector<std::vector<double>> &lattice);

  /**
   * NS distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f;

  /**
   * CDE distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g;

 private:
  // 6  2  5  ^
  //  \ | /   |
  // 3--0--1  |
  //  / | \   |
  // 7  4  8  +------->
  enum Directions {
    E = 1,
    N,
    W,
    S,
    NE,
    NW,
    SW,
    SE
  };

  /**
   * Total duration of the simulation
   */
  double total_time_;

  // by reference, similar to by pointer
  // https://stackoverflow.com/questions/9285627/is-it-possible-to-pass-derived-
  // classes-by-reference-to-a-function-taking-base-cl
  /**
   * Lattice model which contains the number of dimensions, directions, rows,
   * columns, and the discrete velocity vectors
   */
  LatticeModel &lm_;

  /**
   * Collision model for Navier-Stokes equation, contains the equilibrium
   * distribution functions, density and source terms. Performs collision and
   * forcing
   */
  CollisionNS &ns_;

  /**
   * Collision model for Convection-Diffusion equation, contains the equilibrium
   * distribution functions, density and source terms. Performs collision and
   * forcing
   */
  CollisionCD &cd_;

  /**
   * Indicates if Navier-Stokes is implemented in the simulation
   */
  bool is_ns_;

  /**
   * Indicates if Convection-Diffusion is implemented in the simulation
   */
  bool is_cd_;
};
#endif  // LATTICEBOLTZMANN_HPP_
