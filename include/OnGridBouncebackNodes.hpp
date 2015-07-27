#ifndef ON_GRID_BOUNCE_BACK_NODES_HPP_
#define ON_GRID_BOUNCE_BACK_NODES_HPP_
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "StreamModel.hpp"

class OnGridBouncebackNodes: public BoundaryNodes {
 public:
  OnGridBouncebackNodes(bool is_prestream
    , LatticeModel &lm
    , StreamModel &sm
    , CollisionModel &cm);

  ~OnGridBouncebackNodes() = default;

  void AddNode(std::size_t x, std::size_t y);

  void UpdateNodes(std::vector<std::vector<double>> &df);

 protected:
  StreamModel &sm_;
  CollisionModel &cm_;
};

#endif  // ON_GRID_BOUNCE_BACK_NODES_HPP_
