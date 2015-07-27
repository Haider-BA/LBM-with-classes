#include "BouncebackNodes.hpp"
#include <iostream>
#include <vector>
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "Node.hpp"

BouncebackNodes::BouncebackNodes(bool is_prestream
  , LatticeModel &lm
  , CollisionModel *cm)
  : BoundaryNodes(is_prestream, lm),
    cm_ {cm},
    nodes_ {},
    update_nodes {true}
{}

BouncebackNodes::BouncebackNodes(bool is_prestream
  , LatticeModel &lm
  , StreamModel *sm)
  : BoundaryNodes(is_prestream, lm),
    sm_ {sm},
    nodes_ {},
    update_nodes {false}
{}

void BouncebackNodes::AddNode(std::size_t x
  , std::size_t y)
{
  auto nx = lm_.GetNumberOfColumns();
  auto n = y * nx + x;
  nodes_.push_back(Node(x, y, nx));
  // in C++11 nullptr is implicitly cast to boolean false
  // http://stackoverflow.com/questions/11279715/nullptr-and-checking-if-a-
  // pointer-points-to-a-valid-object
  if (cm_) cm_->AddNodeToSkip(n);
  if (sm_) sm_->AddNodeToBounceback(n);
}

void BouncebackNodes::UpdateNodes(std::vector<std::vector<double>> &df)
{
  if (update_nodes) {
    for (auto node : nodes_) {
      auto temp_node = df[node.n];
      df[node.n][E] = temp_node[W];
      df[node.n][N] = temp_node[S];
      df[node.n][W] = temp_node[E];
      df[node.n][S] = temp_node[N];
      df[node.n][NE] = temp_node[SW];
      df[node.n][NW] = temp_node[SE];
      df[node.n][SW] = temp_node[NE];
      df[node.n][SE] = temp_node[NW];
    }  // node
  }
}
