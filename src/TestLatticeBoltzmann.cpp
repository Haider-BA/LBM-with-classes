#include "UnitTest++.h"
#include "Lattice.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <exception>

TEST(DefaultLattice)
{
  Lattice lattice;
  CHECK_EQUAL(2, lattice.GetNumberOfDimensions());
  CHECK_EQUAL(9, lattice.GetNumberOfDiscreteVelocities());
  CHECK_EQUAL(0, lattice.GetNumberOfRows());
  CHECK_EQUAL(0, lattice.GetNumberOfColumns());
}

TEST(NoTypeLattice)
{
  try {
    Lattice lattice(0, 2, 4);
    Lattice lattice2(1, 2, 4);
    Lattice lattice3;
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
    lattice.Print();
    std::cout << "\n";
    lattice2.Print();
    lattice3.Print();
  }
  catch (int e) {
    std::cout << e << "\n";
  }
}
