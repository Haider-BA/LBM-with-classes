#ifndef BOUNDARY_NODES_HPP_
#define BOUNDARY_NODES_HPP_
#include "LatticeModel.hpp"

class BoundaryNodes {
 public:
  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  BoundaryNodes(bool is_prestream
    , bool is_during_stream
    , LatticeModel &lm);

  virtual ~BoundaryNodes() = default;

  virtual void UpdateNodes(std::vector<std::vector<double>> &df
    , bool is_modify_stream) = 0;

  bool prestream;

  bool during_stream;

 protected:
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

  LatticeModel &lm_;
};
#endif // BOUNDARY_NODES_HPP_
