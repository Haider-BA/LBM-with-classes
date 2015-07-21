#ifndef LATTICE_BOLTZMANN_HPP_
#define LATTICE_BOLTZMANN_HPP_
#include <vector>
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "LatticeModel.hpp"
#include "StreamModel.hpp"

class LatticeBoltzmann {
 public:
  /**
   * Constructor: Creates lattice
   * \param t_total total time of simulation
   * \param obstacles_pos lattice containing position of obstacles
   * \param is_ns NS equation toggle
   * \param is_cd CD equation toggle
   * \param is_taylor workaround to toggle bounceback boundary condition for
   *        left and right edge for taylor analytical test
   * \param is_instant instantaneous source togle
   * \param has_obstacles obstacles toggle
   * \param lm lattice model to take care of dimensions, directions, number of
   *        rows, number of columns, space step, time step, omega, discrete
   *        velocity vectors
   * \param ns Collision model for NS equation
   * \param cd Collision model for CD equation
   */
  LatticeBoltzmann(LatticeModel &lm
    , CollisionModel &cm
    , StreamModel &sm);

  ~LatticeBoltzmann() = default;

  LatticeBoltzmann(const LatticeBoltzmann&) = default;

  void AddBoundaryNodes(BoundaryNodes *bn);

  /**
   * Performs one cycle of evolution equation, computes the relevant macroscopic
   * properties such as velocity and density
   */
  void TakeStep();

  /**
   * NS distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> df;

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

  // input parameters
  // LatticeModel to take care of dims, dirs, rows, cols and discrete e vectors
  // by reference, similar to by pointer
  // https://stackoverflow.com/questions/9285627/is-it-possible-to-pass-derived-
  // classes-by-reference-to-a-function-taking-base-cl
  LatticeModel &lm_;

  // Collision models to take care of eq, rho, source
  CollisionModel &cm_;

  StreamModel &sm_;

  // vector of boundary nodes pointer
  std::vector<BoundaryNodes*> bn_;
};
#endif  // LATTICE_BOLTZMANN_HPP_
