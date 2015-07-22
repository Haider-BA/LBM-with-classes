#include "BoundaryNodes.hpp"
#include "LatticeModel.hpp"

BoundaryNodes::BoundaryNodes(bool is_prestream
  , LatticeModel &lm)
  : prestream {is_prestream},
    lm_ (lm)
{}
