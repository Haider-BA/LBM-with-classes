#include "BoundaryNodes.hpp"
#include "LatticeModel.hpp"

BoundaryNodes::BoundaryNodes(LatticeModel &lm
  , bool is_prestream)
  : prestream {is_prestream},
    lm_ (lm)
{}
