#include "OnGridBounceBackNodes.hpp"

OnGridBounceBackNodes::OnGridBounceBackNodes(bool is_prestream
  , LatticeModel &lm
  , StreamModel &sm)
  : BoundaryNodes(is_prestream, lm),
    sm_ (sm)
{}

void OnGridBounceBackNodes::AddNode(std::size_t x
  , std::size_t y)
{
  auto nx = lm_.GetNumberOfColumns();
  auto n = y * nx + x;
  sm_.AddNodeToBounceBack(n);
}

void OnGridBounceBackNodes::UpdateNodes(std::vector<std::vector<double>> &df)
{}
