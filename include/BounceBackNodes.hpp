#ifndef BOUNCE_BACK_NODES_HPP_
#define BOUNCE_BACK_NODES_HPP_
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"

class BounceBackNodes: public BoundaryNodes {
 public:
  /**
   * On-grid bounce-back nodes
   * \param
   * \param
   * \return
   *
   */
  BounceBackNodes(bool is_prestream
    , CollisionModel &cm
    , LatticeModel &lm);

  ~BounceBackNodes() = default;

  void AddNode(std::size_t x, std::size_t y);

  void UpdateNodes(std::vector<std::vector<double>> &df);

 protected:
  CollisionModel &cm_;
};
#endif // BOUNCE_BACK_NODES_HPP_
