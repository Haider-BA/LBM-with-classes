#include "BoundaryNodes.hpp"
#include "LatticeModel.hpp"

BoundaryNodes::BoundaryNodes(bool is_prestream
  , bool is_during_stream
  , LatticeModel &lm)
  : prestream {is_prestream},
    during_stream {is_during_stream},
    lm_ {lm}
{}
