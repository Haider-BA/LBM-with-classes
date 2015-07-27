#ifndef BOUNCE_BACK_NODES_HPP_
#define BOUNCE_BACK_NODES_HPP_
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "Node.hpp"
#include "StreamModel.hpp"

class BouncebackNodes: public BoundaryNodes {
 public:
  /**
   * Creates a full-way bounceback nodes
   * \param
   * \param
   * \return
   */
  BouncebackNodes(bool is_prestream
    , LatticeModel &lm
    , CollisionModel *cm);

  /**
   * Creates a half-way bounceback node
   * \param
   * \param
   */
  BouncebackNodes(bool is_prestream
    , LatticeModel &lm
    , StreamModel *sm);

  /**
   * Override copy constructor due to -Weffc warnings
   */
  BouncebackNodes(const BouncebackNodes&) = default;

  /**
   * Override copy assignment due to -Weffc warnings
   */
  BouncebackNodes& operator= (const BouncebackNodes&) = default;

  ~BouncebackNodes() = default;

  void AddNode(std::size_t x, std::size_t y);

  void UpdateNodes(std::vector<std::vector<double>> &df);

 protected:
  CollisionModel* cm_ = nullptr;

  StreamModel* sm_ = nullptr;

  std::vector<Node> nodes_;

  bool update_nodes;
};
#endif // BOUNCE_BACK_NODES_HPP_
