#include "ZouHePressureNodes.hpp"

ZouHePressureNodes::ZouHePressureNodes(LatticeModel &lm
  , CollisionModel &cm)
  : BoundaryNodes(false, false, lm),
    nodes {},
    cm_ (cm)
{}
