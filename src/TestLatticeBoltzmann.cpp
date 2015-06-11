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
    CHECK_THROW(lattice.Print2D(), std::runtime_error);
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
  u.Init(u0);
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
  std::vector<double> u0 = {1., 2., 3.};
  CHECK_THROW(u.Init(u0), std::runtime_error);
}
/*
TEST(InitUWithWrongLatticeType)
{
  Lattice u(0, 20, 40);
  std::vector<double> u0 = {1., 2.};
  CHECK_THROW(u.Init(u0), std::runtime_error);
}
*/
TEST(CalculateEquilibrium)
{
  double precision = 0.001;
  std::size_t ny = 4;
  std::size_t nx = 3;
  std::size_t nc = 9;
  double density_f = 1.;
  double dx = 1.;
  double dt = 0.001;
  double c = dx / dt;
  std::vector<double> u0 = {1., 1.};
  Lattice f_eq(0, density_f, ny, nx);
  Lattice u(1, ny, nx);
  std::vector<double> ans = {0.444,
      0.111, 0.111, 0.111, 0.111,
      0.028, 0.028, 0.028, 0.028};
  u.Init(u0);
  f_eq.ComputeEq(u, c);
  for (auto y = 0u; y < (nx + 2) * (ny + 2); y += nx + 2) {
    for (auto x = 0u; x < nx + 2; ++x) {
      for (auto i = 0u; i < nc; ++i) {
        if (y == 0 || y == (nx + 2) * (ny + 1) || x == 0 || x == nx + 1) {
          CHECK_CLOSE(0, f_eq.values_[y + x][i], precision);
        }
        else {
          CHECK_CLOSE(ans[i], f_eq.values_[y + x][i], precision);
        }
      }  // i
    }  // x
  }  // y
}

TEST(InitFWithValueFromFEqNode)
{
  double precision = 0.001;
  std::size_t ny = 4;
  std::size_t nx = 3;
  std::size_t nc = 9;
  double density_f = 1.;
  double dx = 1.;
  double dt = 0.001;
  double c = dx / dt;
  std::vector<double> u0 = {123., 321.};
  Lattice f_eq(0, density_f, ny, nx);
  Lattice f(0, ny, nx);
  Lattice u(1, ny, nx);
  u.Init(u0);
  f_eq.ComputeEq(u, c);
  f.Init(f_eq.values_[nx + 3]);
  f.Print2D();
  for (auto y = 0u; y < (nx + 2) * (ny + 2); y += nx + 2) {
    for (auto x = 0u; x < nx + 2; ++x) {
      for (auto i = 0u; i < nc; ++i) {
        CHECK_CLOSE(f_eq.values_[y + x][i], f.values_[y + x][i], precision);
      }  // i
    }  // x
  }  // y
}
