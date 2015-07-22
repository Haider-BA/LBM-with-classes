#include "BounceBackNodes.hpp"
#include <iostream>
#include <vector>
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "Node.hpp"

BounceBackNodes::BounceBackNodes(bool is_prestream
  , CollisionModel &cm
  , LatticeModel &lm)
  : BoundaryNodes(is_prestream, lm),
    cm_ (cm),
    nodes_ {}
{}

void BounceBackNodes::AddNode(std::size_t x
  , std::size_t y)
{
  auto nx = lm_.GetNumberOfColumns();
  auto n = y * nx + x;
  nodes_.push_back(Node(x, y, nx));
  cm_.AddNodeToSkip(n);
}

void BounceBackNodes::UpdateNodes(std::vector<std::vector<double>> &df)
{
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
  }
}
