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
   * Constructor: Creates a LatticeBoltzmann object
   * \param lm lattice model which contains information on the number of rows,
   *        columns, dimensions, discrete directions and lattice velocity
   * \param cm collision model used by the lattice: Convection-diffusion,
   *        Navier-Stokes and Navier-Stokes with force
   * \param sm stream mode used by the lattice: Periodic stream, non-periodic
   *        streaming
   */
  LatticeBoltzmann(LatticeModel &lm
    , CollisionModel &cm
    , StreamModel &sm);

  /**
   * Destructor
   */
  ~LatticeBoltzmann() = default;

  /**
   * Copy constructor
   */
  LatticeBoltzmann(const LatticeBoltzmann&) = default;

  /**
   * Adds a boundary condition to the lattice
   * \param bn pointer to the boundary condition to be added
   */
  void AddBoundaryNodes(BoundaryNodes *bn);

  /**
   * Performs one cycle of evolution equation, computes the relevant macroscopic
   * properties such as velocity and density
   */
  void TakeStep();

  /**
   * Lattice distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> df;

 private:
  // 6  2  5  ^
  //  \ | /   |
  // 3--0--1  |
  //  / | \   |
  // 7  4  8  +------->
  /**
   * Enumeration of discrete directions to be used with lattice distribution
   * functions
   */
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

  // by reference, similar to by pointer
  // https://stackoverflow.com/questions/9285627/is-it-possible-to-pass-derived-
  // classes-by-reference-to-a-function-taking-base-cl
  /**
   * Lattice model which contains information on the number of rows, columns,
   * dimensions, discrete directions and lattice velocity
   */
  LatticeModel &lm_;

  /**
   * Collision models to perform the collision step based on the model chosen
   */
  CollisionModel &cm_;

  /**
   * Stream model to perform the streaming step based on the model chosen
   */
  StreamModel &sm_;

  /**
   * Pointers to boundary conditions in the lattice stored in a vector.
   * References cannot be used as it is not possible to store a vector of
   * references
   */
  std::vector<BoundaryNodes*> bn_;
};
#endif  // LATTICE_BOLTZMANN_HPP_
