#ifndef BOUNDARY_NODES_HPP_
#define BOUNDARY_NODES_HPP_
#include "LatticeModel.hpp"

class BoundaryNodes {
 public:
  /**
   * Creates the base for all boundary nodes, contains boolean toggles for
   * prestream and during stream functions so that a single function call can be
   * used in the LatticeBoltzmann TakeStep() method
   * \param is_prestream Boolean toggle to indicate if boundary condition has
   *        any functions to be executed prestream
   * \param is_during_stream Boolean toggle to indicate if boundary condition
   *        has any functions to be executed during stream, or after the
   *        Stream() method based on implementation in LatticeBoltzmann class
   * \param lm reference to LatticeModel which contains information on the
   *        number of rows, columns, dimensions, discrete directions and lattice
   *        velocity
   */
  BoundaryNodes(bool is_prestream
    , bool is_during_stream
    , LatticeModel &lm);

  /**
   * Virtual destructor since we are deriving from this class
   */
  virtual ~BoundaryNodes() = default;

  /**
   * Pure virtual function for the boundary conditions to implement on how the
   * boundary nodes are updated
   * \param df lattice distribution functions stored row-wise in a 2D vector
   * \param is_modify_stream Boolean toggle for half-way bounceback nodes to
   *        perform functions after streaming
   */
  virtual void UpdateNodes(std::vector<std::vector<double>> &df
    , bool is_modify_stream) = 0;

  /**
   * Boolean toggle to indicate if boundary condition occurs before streaming
   */
  bool prestream;

  /**
   * Boolean toggle to indicate if boundary condition occurs during stream, or
   * after Stream() function based on implementation
   */
  bool during_stream;

 protected:
  /**
   * Enumeration for discrete directions to be used with df
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

  /**
   * Reference to LatticeModel which contains information on number of rows,
   * columns, dimensions, discrete directions and lattice velocity
   */
  LatticeModel &lm_;
};
#endif // BOUNDARY_NODES_HPP_
