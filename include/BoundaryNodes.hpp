#ifndef BOUNDARY_NODES_HPP_
#define BOUNDARY_NODES_HPP_
#include "LatticeModel.hpp"

class BoundaryNodes {
 public:
  BoundaryNodes(LatticeModel &lm
    , bool is_prestream);

  virtual ~BoundaryNodes() = default;

  virtual void UpdateNodes(std::vector<std::vector<double>> &df) = 0;

  bool prestream;

 protected:
  LatticeModel &lm_;
};
#endif // BOUNDARY_NODES_HPP_
