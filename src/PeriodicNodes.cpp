#include "PeriodicNodes.hpp"
#include <vector>
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"

PeriodicNodes::PeriodicNodes(bool is_prestream
  , CollisionModel &cm
  , LatticeModel &lm)
  : BoundaryNodes(is_prestream, lm),
    cm_ (cm)
{}

void PeriodicNodes::AddNode(std::size_t x
  , std::size_t y)
{

}
