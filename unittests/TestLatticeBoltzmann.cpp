#include <stdexcept>  // runtime_error
#include "LatticeBoltzmann.hpp"
#include "UnitTest++.h"

SUITE(TestException)
{
TEST(ConstructiorException)
{
  CHECK_THROW(LatticeBoltzmann lbm, std::runtime_error);
}
}
