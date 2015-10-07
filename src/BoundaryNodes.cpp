#include "BoundaryNodes.hpp"
#include <vector>
#include "LatticeModel.hpp"

// have to use parenthesis for reference initializing in initializer list due to
// bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=50025
BoundaryNodes::BoundaryNodes(bool is_prestream
  , bool is_during_stream
  , LatticeModel &lm)
  : prestream {is_prestream},
    during_stream {is_during_stream},
    position {},
    lm_ (lm)
{}
