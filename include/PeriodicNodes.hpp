#ifndef PERIODIC_NODES_HPP_
#define PERIODIC_NODES_HPP_
#include "BoundaryNodes.hpp"

class PeriodicNodes: public BoundaryNodes {
 public:
  PeriodicNodes(bool is_prestream
    , CollisionModel &cm
    , LatticeModel &lm);

  ~PeriodicNodes() = default;

  void AddNode(std::size_t x, std::size_t y);

  void UpdateNodes(std::vector<std::vector<double>> &df);

 protected:
  CollisionModel cm_;
};
#endif // PERIODIC_NODES_HPP_
