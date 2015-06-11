#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Lattice.hpp"
#include "UnitTest++.h"

TEST(DefaultLattice)
{
  Lattice lattice;
  CHECK_EQUAL(2, lattice.GetNumberOfDimensions());
  CHECK_EQUAL(9, lattice.GetNumberOfDiscreteVelocities());
  CHECK_EQUAL(0, lattice.GetNumberOfRows());
  CHECK_EQUAL(0, lattice.GetNumberOfColumns());
}

TEST(LatticeWithIncorrectType)
{
  CHECK_THROW(Lattice lattice(-1,2,4), std::runtime_error);
}

TEST(PrintIncorrectLatticeType)
{

    Lattice lattice;
    CHECK_THROW(lattice.Print(), std::runtime_error);
//    Lattice lattice2(2, 2, 4);
//    CHECK_THROW(lattice2.Print(), std::runtime_error);
}

TEST(Depth9Lattice)
{
  Lattice lattice(0, 2, 4);
  std::size_t nx = lattice.GetNumberOfColumns();
  std::size_t ny = lattice.GetNumberOfRows();
  CHECK_EQUAL(2, lattice.GetNumberOfDimensions());
  CHECK_EQUAL(9, lattice.GetNumberOfDiscreteVelocities());
  CHECK_EQUAL(2, ny);
  CHECK_EQUAL(4, nx);
  for (auto y = 0u; y < (nx + 2) * (ny + 2); y += nx + 2) {
    for (auto x = 0u; x < nx + 2; ++x) {
      CHECK_EQUAL(9, lattice.values_[y + x].size());
    } // x
  } // y
}

TEST(Depth2Lattice)
{
  Lattice lattice(1, 2, 4);
  std::size_t nx = lattice.GetNumberOfColumns();
  std::size_t ny = lattice.GetNumberOfRows();
  CHECK_EQUAL(2, lattice.GetNumberOfDimensions());
  CHECK_EQUAL(9, lattice.GetNumberOfDiscreteVelocities());
  CHECK_EQUAL(2, ny);
  CHECK_EQUAL(4, nx);
  for (auto y = 0u; y < (nx + 2) * (ny + 2); y += nx + 2) {
    for (auto x = 0u; x < nx + 2; ++x) {
      CHECK_EQUAL(2, lattice.values_[y + x].size());
    } // x
  } // y
}

TEST(InitUWithProperValues)
{
  Lattice u(1, 20, 40);
  std::size_t nx = u.GetNumberOfColumns();
  std::size_t ny = u.GetNumberOfRows();
  std::size_t nd = u.GetNumberOfDimensions();
  std::vector<double> u0 = {1., 2.};
  u.InitU(u0);
  for (auto y = nx + 2; y < (nx + 2) * (ny + 1); y += nx + 2) {
    for (auto x = 1u; x < nx + 1; ++x) {
      for (auto d = 0u; d < nd; ++d) {
        CHECK_EQUAL(u0[d], u.values_[y + x][d]);
      }
    } // x
  } // y
}

TEST(InitUWithWrongValues)
{
  Lattice u(1, 20, 40);
  std::size_t nx = u.GetNumberOfColumns();
  std::size_t ny = u.GetNumberOfRows();
  std::size_t nd = u.GetNumberOfDimensions();
  std::vector<double> u0 = {1., 2., 3.};
  CHECK_THROW(u.InitU(u0), std::runtime_error);
}
