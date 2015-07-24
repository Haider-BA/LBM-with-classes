#ifndef ON_GRID_BOUNCE_BACK_NODES_HPP_
#define ON_GRID_BOUNCE_BACK_NODES_HPP_
#include "BoundaryNodes.hpp"
#include "StreamModel.hpp"

class OnGridBounceBackNodes: public BoundaryNodes {
 public:
  OnGridBounceBackNodes(bool is_prestream
    , LatticeModel &lm
    , StreamModel &sm);

  ~OnGridBounceBackNodes() = default;

  void AddNode(std::size_t x, std::size_t y);

  void UpdateNodes(std::vector<std::vector<double>> &df);

 protected:
  StreamModel &sm_;
};

#endif  // ON_GRID_BOUNCE_BACK_NODES_HPP_
