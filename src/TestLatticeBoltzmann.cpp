#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Lattice.hpp"
#include "UnitTest++.h"

enum DiscreteDirections {
  E = 1,
  N,
  W,
  S,
  NE,
  NW,
  SW,
  SE,
};

TEST(DefaultLattice)
{
  Lattice lattice;
  CHECK_EQUAL(0u, lattice.GetNumberOfDimensions());
  CHECK_EQUAL(0u, lattice.GetNumberOfDiscreteVelocities());
  CHECK_EQUAL(0u, lattice.GetNumberOfRows());
  CHECK_EQUAL(0u, lattice.GetNumberOfColumns());
}

TEST(ZeroValuesInLatticeDeclaration)
{
  CHECK_THROW(Lattice lattice(0, 9, 2, 4, 1, 1, true, true),
      std::runtime_error);
  CHECK_THROW(Lattice lattice(2, 0, 2, 4, 1, 1, true, true),
      std::runtime_error);
  CHECK_THROW(Lattice lattice(2, 9, 0, 4, 1, 1, true, true),
      std::runtime_error);
  CHECK_THROW(Lattice lattice(2, 9, 2, 0, 1, 1, true, true),
      std::runtime_error);
  CHECK_THROW(Lattice lattice(2, 9, 2, 4, 0, 1, true, true),
      std::runtime_error);
  CHECK_THROW(Lattice lattice(2, 9, 2, 4, 1, 0, true, true),
      std::runtime_error);
  CHECK_THROW(Lattice lattice(2, 9, 2, 4, 1, 1, false, false),
      std::runtime_error);
}

TEST(NormalLatticeBeforeInit)
{
  Lattice lattice(2, 9, 2, 4, 1, 1, true, true);
  std::size_t ny = lattice.GetNumberOfRows();
  std::size_t nx = lattice.GetNumberOfColumns();
  std::size_t nd = lattice.GetNumberOfDimensions();
  std::size_t nc = lattice.GetNumberOfDiscreteVelocities();
  std::size_t lattice_size = ny * nx;
  CHECK_EQUAL(2u, nd);
  CHECK_EQUAL(9u, nc);
  CHECK_EQUAL(2u, ny);
  CHECK_EQUAL(4u, nx);
  for (auto n = 0u; n < lattice_size; ++n) {
    CHECK_EQUAL(nc, lattice.f_[n].size());
    CHECK_EQUAL(nc, lattice.g_[n].size());
    CHECK_EQUAL(nc, lattice.f_eq_[n].size());
    CHECK_EQUAL(nc, lattice.g_eq_[n].size());
    CHECK_EQUAL(nd, lattice.src_f_[n].size());
    CHECK_EQUAL(nd, lattice.u_[n].size());
    for (auto i = 0u; i < nc; ++i) {
      CHECK_EQUAL(0.0, lattice.f_[n][i]);
      CHECK_EQUAL(0.0, lattice.g_[n][i]);
      CHECK_EQUAL(0.0, lattice.f_eq_[n][i]);
      CHECK_EQUAL(0.0, lattice.g_eq_[n][i]);
    }  // i
    for (auto d = 0u; d < nd; ++d) {
      CHECK_EQUAL(0.0, lattice.src_f_[n][d]);
      CHECK_EQUAL(0.0, lattice.u_[n][d]);
    }  // d
    CHECK_EQUAL(0.0, lattice.src_g_[n]);
    CHECK_EQUAL(0.0, lattice.rho_f_[n]);
    CHECK_EQUAL(0.0, lattice.rho_g_[n]);
  }  // n
}

TEST(Init1DLattices)
{
  Lattice lattice(2, 9, 2, 4, 1, 1, true, true);
  const double initial_rho = 1.0;
  lattice.Init(lattice.rho_f_, initial_rho);
  lattice.Init(lattice.rho_g_, initial_rho);
  for (auto lat : lattice.rho_f_) {
    CHECK_EQUAL(initial_rho, lat);
  }
  for (auto lat : lattice.rho_g_) {
    CHECK_EQUAL(initial_rho, lat);
  }
}

TEST(Init2DLatticesWithValues)
{
  Lattice lattice(2, 9, 2, 4, 1, 1, true, true);
  std::vector<double> initial_u = {1.0, 2.0};
  lattice.Init(lattice.u_, initial_u);
  for (auto lat : lattice.u_) {
    CHECK_EQUAL(initial_u[0], lat[0]);
    CHECK_EQUAL(initial_u[1], lat[1]);
  }
}

TEST(Init2DLatticesWithValuesWrongDimensions)
{
  Lattice lattice(2, 9, 2, 4, 1, 1, true, true);
  std::vector<double> initial_u = {1.0, 2.0, 3.0};
  CHECK_THROW(lattice.Init(lattice.u_, initial_u), std::runtime_error);
}

TEST(Init2DLatticesWithLattice)
{
  Lattice lattice(2, 9, 2, 4, 1, 1, true, true);
  auto ny = lattice.GetNumberOfRows();
  auto nx = lattice.GetNumberOfColumns();
  std::vector<double> node = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<std::vector<double>> initial_g(nx * ny, node);
  lattice.Init(lattice.g_, initial_g);
  for (auto lat : lattice.g_) {
    int index = 0;
    for (auto i : lat) CHECK_EQUAL(node[index++], i);
  }  // lat
}

TEST(Init2DLatticesWithLatticeWrongDimensions)
{
  Lattice lattice(2, 9, 2, 4, 1, 1, true, true);
  auto ny = lattice.GetNumberOfRows();
  auto nx = lattice.GetNumberOfColumns();
  std::vector<double> node = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<std::vector<double>> initial_g(nx * (ny + 2), node);
  CHECK_THROW(lattice.Init(lattice.g_, initial_g), std::runtime_error);
}

TEST(Init2DLatticesWithLatticeWrongDepth)
{
  Lattice lattice(2, 9, 2, 4, 1, 1, true, true);
  auto ny = lattice.GetNumberOfRows();
  auto nx = lattice.GetNumberOfColumns();
  auto nc = lattice.GetNumberOfDiscreteVelocities();
  for (auto n = 0u; n < 100; ++n) {
    if (n == nc) continue;
    std::vector<double> node(n, 0.0);
    std::vector<std::vector<double>> initial_g((nx + 2) * (ny + 2), node);
    CHECK_THROW(lattice.Init(lattice.g_, initial_g), std::runtime_error);
  }
}

TEST(CalculateEquilibrium)
{
  double precision = 0.001;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 4;
  std::size_t nx = 3;
  double density_g = 1.;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  std::vector<double> u0 = {1., 1.};
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  std::vector<double> ans = {0.444,
      0.111, 0.111, 0.111, 0.111,
      0.028, 0.028, 0.028, 0.028};
  lattice.Init(lattice.u_, u0);
  lattice.Init(lattice.rho_g_, density_g);
  lattice.ComputeEq(lattice.g_eq_, lattice.rho_g_);

  for (auto lat_eq : lattice.g_eq_) {
    for (auto i = 0u; i < nc; ++i) {
      CHECK_CLOSE(ans[i], lat_eq[i], precision);
    }  // i
  }  // lat_eq
}

TEST(InitFWithValueFromFEqNode)
{
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 4;
  std::size_t nx = 3;
  double density_g = 1.;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  std::vector<double> u0 = {123., 321.};
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  lattice.Init(lattice.u_, u0);
  lattice.Init(lattice.rho_g_, density_g);
  lattice.ComputeEq(lattice.g_eq_, lattice.rho_g_);
  lattice.Init(lattice.g_, lattice.g_eq_);
  for (auto n = 0u; n < nx * ny; ++n) {
    for (auto i = 0u; i < nc; ++i) {
      CHECK_EQUAL(lattice.g_eq_[n][i], lattice.g_[n][i]);
    }  // i
  }  // n
}

TEST(InitSingleSource)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 4;
  std::size_t nx = 5;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  unsigned x_pos = 4;
  unsigned y_pos = 2;
  std::vector<std::vector<unsigned>> src_position = {{x_pos, y_pos}};
  std::vector<double> src_strength = {2};
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  lattice.InitSrc(lattice.src_g_, src_position, src_strength);
  for (auto n = 0u; n < nx * ny; ++n) {
    if (n == y_pos * nx + x_pos) {
      CHECK_CLOSE(src_strength[0], lattice.src_g_[n], zero_tol);
    }
    else {
      CHECK_CLOSE(0.0, lattice.src_g_[n], zero_tol);
    }
  }  // n
}

TEST(InitMultipleSource)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 10;
  std::size_t nx = 10;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  unsigned x_pos = 4;
  unsigned y_pos = 2;
  std::vector<std::vector<unsigned>> src_position = {{x_pos, y_pos},
      {x_pos + 1, y_pos + 2}};
  std::vector<double> src_strength = {2, 3.4};
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  lattice.InitSrc(lattice.src_g_, src_position, src_strength);
  for (auto n = 0u; n < nx * ny; ++n) {
    if (n == y_pos * nx + x_pos) {
      CHECK_CLOSE(src_strength[0], lattice.src_g_[n], zero_tol);
    }
    else if (n == (y_pos + 2) * nx + x_pos + 1) {
      CHECK_CLOSE(src_strength[1], lattice.src_g_[n], zero_tol);
    }
    else {
      CHECK_CLOSE(0.0, lattice.src_g_[n], zero_tol);
    }
  }  // n
}

TEST(InitSingle2DSource)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 4;
  std::size_t nx = 5;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  unsigned x_pos = 4;
  unsigned y_pos = 2;
  std::vector<std::vector<unsigned>> src_position = {{x_pos, y_pos}};
  std::vector<std::vector<double>> src_strength = {{2, 3}};
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  lattice.InitSrc(lattice.src_f_, src_position, src_strength);
  for (auto n = 0u; n < nx * ny; ++n) {
    if (n == y_pos * nx + x_pos) {
      CHECK_CLOSE(src_strength[0][0], lattice.src_f_[n][0], zero_tol);
      CHECK_CLOSE(src_strength[0][1], lattice.src_f_[n][1], zero_tol);
    }
    else {
      CHECK_CLOSE(0.0, lattice.src_f_[n][0], zero_tol);
      CHECK_CLOSE(0.0, lattice.src_f_[n][1], zero_tol);
    }
  }  // n
}

TEST(InitMultiple2DSource)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 10;
  std::size_t nx = 11;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  unsigned x_pos = 4;
  unsigned y_pos = 2;
  std::vector<std::vector<unsigned>> src_position = {{x_pos, y_pos},
      {x_pos + 2, y_pos + 1}};
  std::vector<std::vector<double>> src_strength = {{2, 3}, {2.3, 4.4}};
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  lattice.InitSrc(lattice.src_f_, src_position, src_strength);
  for (auto n = 0u; n < nx * ny; ++n) {
    if (n == y_pos * nx + x_pos) {
      CHECK_CLOSE(src_strength[0][0], lattice.src_f_[n][0], zero_tol);
      CHECK_CLOSE(src_strength[0][1], lattice.src_f_[n][1], zero_tol);
    }
    else if (n == (y_pos + 1) * nx + x_pos + 2) {
      CHECK_CLOSE(src_strength[1][0], lattice.src_f_[n][0], zero_tol);
      CHECK_CLOSE(src_strength[1][1], lattice.src_f_[n][1], zero_tol);
    }
    else {
      CHECK_CLOSE(0.0, lattice.src_f_[n][0], zero_tol);
      CHECK_CLOSE(0.0, lattice.src_f_[n][1], zero_tol);
    }
  }  // n
}

// TODO CHECK_THROW for InitSrc

TEST(BoundaryConditionPeriodic)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 3;
  std::size_t nx = 4;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  std::vector<double> nums = {0, 1, 0, 2, 4, 5, 6, 7, 8};
  lattice.Init(lattice.g_, nums);
  lattice.BoundaryCondition(lattice.g_, lattice.boundary_g_);
  for (auto y = 0u; y < ny; ++y) {
    auto n = y * (nx);
    CHECK_CLOSE(lattice.g_[n + nx - 1][E], lattice.boundary_g_[y][E], zero_tol);
    CHECK_CLOSE(lattice.g_[n + nx - 1][NE], lattice.boundary_g_[y][NE],
        zero_tol);
    CHECK_CLOSE(lattice.g_[n + nx - 1][SE], lattice.boundary_g_[y][SE],
        zero_tol);
    CHECK_CLOSE(lattice.g_[n][W], lattice.boundary_g_[y + ny][W], zero_tol);
    CHECK_CLOSE(lattice.g_[n][NW], lattice.boundary_g_[y + ny][NW], zero_tol);
    CHECK_CLOSE(lattice.g_[n][SW], lattice.boundary_g_[y + ny][SW], zero_tol);
  }  // y
}

TEST(BoundaryConditionBounceBack)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 3;
  std::size_t nx = 4;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  std::vector<double> nums = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  double value = 0;
  for (auto &lat : lattice.g_) lat = std::vector<double>(nc, value++);
  lattice.BoundaryCondition(lattice.g_, lattice.boundary_g_);
  auto top = 2 * ny;
  auto bottom = 2 * ny + nx;
  for (auto x = 0u; x < nx; ++x) {
    auto n = (ny - 1) * (nx);
    CHECK_CLOSE(lattice.g_[x + n][N], lattice.boundary_g_[top + x][S],
        zero_tol);
    CHECK_CLOSE(lattice.g_[x][S], lattice.boundary_g_[bottom + x][N],
        zero_tol);
    if (x == 0) {
      CHECK_CLOSE(lattice.g_[n + ny][NE], lattice.boundary_g_[top + x][SW],
          zero_tol);
      CHECK_CLOSE(lattice.g_[nx - 1][SE], lattice.boundary_g_[bottom + x][NW],
          zero_tol);
    }
    else {
      CHECK_CLOSE(lattice.g_[x + n - 1][NE], lattice.boundary_g_[top + x][SW],
          zero_tol);
      CHECK_CLOSE(lattice.g_[x - 1][SE], lattice.boundary_g_[bottom + x][NW],
          zero_tol);
    }
    if (x == nx - 1) {
      CHECK_CLOSE(lattice.g_[n][NW], lattice.boundary_g_[top + x][SE],
          zero_tol);
      CHECK_CLOSE(lattice.g_[0][SW], lattice.boundary_g_[bottom + x][NE],
          zero_tol);
    }
    else {
      CHECK_CLOSE(lattice.g_[x + n + 1][NW], lattice.boundary_g_[top + x][SE],
          zero_tol);
      CHECK_CLOSE(lattice.g_[x + 1][SW], lattice.boundary_g_[bottom + x][NE],
          zero_tol);
    }
  }
}

TEST(BoundaryConditionCorner)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 5;
  std::size_t nx = 4;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  std::vector<double> nums = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  lattice.Init(lattice.g_, nums);
  lattice.BoundaryCondition(lattice.g_, lattice.boundary_g_);
  auto corner = 2 * nx + 2 * ny;
  CHECK_CLOSE(lattice.g_[0][SW], lattice.boundary_g_[corner][NE], zero_tol);
  CHECK_CLOSE(lattice.g_[nx - 1][SE], lattice.boundary_g_[corner + 1][NW],
      zero_tol);
  CHECK_CLOSE(lattice.g_[(ny - 1) * nx][NW],
      lattice.boundary_g_[corner + 2][SE], zero_tol);
  CHECK_CLOSE(lattice.g_[ny * nx - 1][NE], lattice.boundary_g_[corner + 3][SW],
      zero_tol);
}

TEST(StreamHorizontal)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 6;
  std::size_t nx = 4;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  std::vector<double> first_three = {1, 2, 1, 0, 1, 2, 0, 0, 2};
  std::vector<double> second_three = {4, 5, 4, 3, 4, 5, 3, 3, 5};
  std::vector<std::vector<double>> nodes = {first_three, second_three};
  std::vector<double> first_result = {1, 5, 1, 3, 1, 5, 3, 3, 5};
  std::vector<double> second_result = {4, 2, 4, 0, 4, 2, 0, 0, 2};
  int counter = 0;
  for (auto &lat : lattice.g_) lat = nodes[counter++ % 2];
  lattice.BoundaryCondition(lattice.g_, lattice.boundary_g_);
  for (auto n = 0u; n < lattice.boundary_g_.size(); ++n) {
    if (n < ny) {
      lattice.boundary_g_[n] = nodes[1];
    }
    else if (n < 2 * ny) {
      lattice.boundary_g_[n] = nodes[ny % 2];
    }
    else if (n < 2 * ny + 2 * nx) {
      lattice.boundary_g_[n] = nodes[n % 2];
    }
    else {
      lattice.boundary_g_[n] = nodes[(n + 1) % 2];
    }
  }  // n
  lattice.g_ = lattice.Stream(lattice.g_, lattice.boundary_g_);
  for (auto lat : lattice.g_) {
    for (auto i = 0u; i < nc; ++i) {
      if ((lat[0] - 1) < zero_tol) {
      CHECK_CLOSE(first_result[i], lat[i], zero_tol);
      }
      else {
        CHECK_CLOSE(second_result[i], lat[i], zero_tol);
      }
    }  // i
  }  // lat
}

TEST(StreamVertical)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 5;
  std::size_t nx = 4;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  std::vector<double> first_three = {1, 1, 0, 1, 2, 0, 0, 2, 2};
  std::vector<double> second_three = {4, 4, 3, 4, 5, 3, 3, 5, 5};
  std::vector<std::vector<double>> nodes = {first_three, second_three};
  std::vector<double> first_result = {1, 1, 3, 1, 5, 3, 3, 5, 5};
  std::vector<double> second_result = {4, 4, 0, 4, 2, 0, 0, 2, 2};
  int counter = 0;
  for (auto &lat : lattice.g_) lat = nodes[(counter++ / nx) % 2];
  lattice.BoundaryCondition(lattice.g_, lattice.boundary_g_);
  for (auto n = 0u; n < lattice.boundary_g_.size(); ++n) {
    if (n < ny) {
      lattice.boundary_g_[n] = nodes[n % 2];
    }
    else if (n < 2 * ny) {
      lattice.boundary_g_[n] = nodes[(n + 1) % 2];
    }
    else if (n < 2 * ny + 2 * nx) {
      lattice.boundary_g_[n] = nodes[1];
    }
    else {
      lattice.boundary_g_[n] = nodes[(n % 2 + (n + 1) % 2) % 2];
    }
  }  // n
  lattice.g_ = lattice.Stream(lattice.g_, lattice.boundary_g_);
  for (auto lat : lattice.g_) {
    for (auto i = 0u; i < nc; ++i) {
      if ((lat[0] - 1) < zero_tol) {
      CHECK_CLOSE(first_result[i], lat[i], zero_tol);
      }
      else {
        CHECK_CLOSE(second_result[i], lat[i], zero_tol);
      }
    }  // i
  }  // lat
}

TEST(StreamDiagonalNESW)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 5;
  std::size_t nx = 4;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  std::vector<double> ones(9, 1);
  std::vector<double> twos(9, 2);
  std::vector<double> threes(9, 3);
  std::vector<std::vector<double>> nodes = {ones, twos, threes};
  std::vector<double> result = {1, 2, 3};
  for (auto y = 0u; y < ny; ++y) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      lattice.g_[n] = nodes[(x % 3 + y % 3) % 3];
    }  // x
  }  // y
  lattice.BoundaryCondition(lattice.g_, lattice.boundary_g_);
  for (auto n = 0u; n < lattice.boundary_g_.size(); ++n) {
    if (n < ny) {
      lattice.boundary_g_[n] = nodes[(n + 2) % 3];
    }
    else if (n < 2 * ny) {
      lattice.boundary_g_[n] = nodes[(n + 2) % 3];
    }
    else if (n < 2 * ny + nx) {
      lattice.boundary_g_[n] = nodes[(n + 1) % 3];
    }
    else if (n < 2 * ny + 2 * nx) {
      lattice.boundary_g_[n] = nodes[n % 3];
    }
    else {
      lattice.boundary_g_[n] = nodes[(n + 1) % 2];
    }
  }  // n
  lattice.g_ = lattice.Stream(lattice.g_, lattice.boundary_g_);
  for (auto y = 0u; y < ny; ++y) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      CHECK_CLOSE(result[(x % 3 + y % 3 + 1) % 3], lattice.g_[n][NE], zero_tol);
      CHECK_CLOSE(result[(x % 3 + y % 3 + 2) % 3], lattice.g_[n][SW], zero_tol);
    }  // x
  }  // y
}

TEST(StreamDiagonalNWSE)
{
  double zero_tol = 1e-20;
  std::size_t nd = 2;
  std::size_t nc = 9;
  std::size_t ny = 5;
  std::size_t nx = 4;
  double dx = 1.;
  double dt = 0.001;
  bool is_cd = true;
  bool is_ns = true;
  Lattice lattice(nd, nc, ny, nx, dx, dt, is_cd, is_ns);
  std::vector<double> ones(9, 1);
  std::vector<double> twos(9, 2);
  std::vector<double> threes(9, 3);
  std::vector<std::vector<double>> nodes = {ones, twos, threes};
  std::vector<double> result = {1, 2, 3};
  for (auto y = 0u; y < ny; ++y) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      lattice.g_[n] = nodes[(2 - x % 3 + y % 3) % 3];
    }  // x
  }  // y
  lattice.BoundaryCondition(lattice.g_, lattice.boundary_g_);
  for (auto n = 0u; n < lattice.boundary_g_.size(); ++n) {
    if (n < ny) {
      lattice.boundary_g_[n] = nodes[n % 3];
    }
    else if (n < 2 * ny) {
      lattice.boundary_g_[n] = nodes[(n + 2) % 3];
    }
    else if (n < 2 * ny + nx) {
      lattice.boundary_g_[n] = nodes[(1 - n) % 3];
    }
    else if (n < 2 * ny + 2 * nx) {
      lattice.boundary_g_[n] = nodes[(2 - n) % 3];
    }
    else {
      lattice.boundary_g_[n] = nodes[2 * ((n + 1) % 2)];
    }
  }  // n
  lattice.g_ = lattice.Stream(lattice.g_, lattice.boundary_g_);
  for (auto y = 0u; y < ny; ++y) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      CHECK_CLOSE(result[(3 - x % 3 + y % 3) % 3], lattice.g_[n][NW], zero_tol);
      CHECK_CLOSE(result[(4 - x % 3 + y % 3) % 3], lattice.g_[n][SE], zero_tol);
    }  // x
  }  // y
}
