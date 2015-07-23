#ifndef ZOU_HE_NODES_HPP_
#define ZOU_HE_NODES_HPP_
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "ValueNode.hpp"

class ZouHeNodes: public BoundaryNodes {
 public:
  ZouHeNodes(bool is_prestream
    , CollisionModel &cm
    , LatticeModel &lm);

  ~ZouHeNodes() = default;

  void AddNode(std::size_t x
    , std::size_t y
    , double u_lid
    , double v_lid);

  void UpdateNodes(std::vector<std::vector<double>> &df);

  void UpdateSide(std::vector<std::vector<double>> &df
    , ValueNode &node);

  void UpdateCorner(std::vector<std::vector<double>> &df
    , std::size_t n);

 protected:
  std::vector<ValueNode> nodes_;
//  std::vector<bool> is_corner_;
//  std::vector<bool> knowns_;

};

#endif  // ZOU_HE_NODES_HPP_
