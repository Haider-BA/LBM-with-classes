#include "StreamModel.hpp"
#include "LatticeModel.hpp"

StreamModel::StreamModel(LatticeModel &lm)
  : lm_ (lm),
    bounce_back {}
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  bounce_back.assign(nx * ny, false);
}

void StreamModel::AddNodeToBounceback(std::size_t n)
{
  bounce_back[n] = true;
}
