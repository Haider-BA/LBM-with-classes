#include "BounceBackNodes.hpp"
#include <vector>
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"

BounceBackNodes::BounceBackNodes(bool is_prestream
  , CollisionModel &cm
  , LatticeModel &lm)
  : BoundaryNodes(is_prestream, lm),
    cm_ (cm)
{}

void BounceBackNodes::AddNode(std::size_t x
  , std::size_t y)
{
  auto nx = lm_.GetNumberOfColumns();
  auto n = y * nx + x;
  coordinates_.push_back(n);
  cm_.AddNodeToSkip(n);
}

void BounceBackNodes::UpdateNodes(std::vector<std::vector<double>> &df)
{
  for (auto coord : coordinates_) {
    auto temp_node = df[coord];
    df[coord][E] = temp_node[W];
    df[coord][N] = temp_node[S];
    df[coord][W] = temp_node[E];
    df[coord][S] = temp_node[N];
    df[coord][NE] = temp_node[SW];
    df[coord][NW] = temp_node[SE];
    df[coord][SW] = temp_node[NE];
    df[coord][SE] = temp_node[NW];
  }
}
