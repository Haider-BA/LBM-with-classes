#include "LatticeBoltzmann.hpp"
#include <stdexcept>  // runtime_error
#include "LatticeModel.hpp"

LatticeBoltzmann::LatticeBoltzmann()
{
  throw std::runtime_error("Please use the other constructor");
}

LatticeBoltzmann::LatticeBoltzmann(const LatticeModel &lm)
{
  throw std::runtime_error("Please use the other constructor");
}
