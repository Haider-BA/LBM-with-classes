#include "OnGridBouncebackNodes.hpp"

OnGridBouncebackNodes::OnGridBouncebackNodes(bool is_prestream
  , LatticeModel &lm
  , StreamModel &sm
  , CollisionModel &cm)
  : BoundaryNodes(is_prestream, lm),
    sm_ (sm),
    cm_ (cm)
{}

void OnGridBouncebackNodes::AddNode(std::size_t x
  , std::size_t y)
{
  auto nx = lm_.GetNumberOfColumns();
  auto n = y * nx + x;
  sm_.AddNodeToBounceback(n);
//  cm_.AddNodeToSkip(n);
}

void OnGridBouncebackNodes::UpdateNodes(std::vector<std::vector<double>> &df)
{}
