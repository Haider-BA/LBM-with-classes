#include "LatticeBoltzmann.hpp"
#include <cmath>  // std::fmod
#include <iomanip>  // std::setprecision
#include <iostream>
#include <stdexcept>  // std::runtime_error
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeModel.hpp"
#include "Printing.hpp"
#include "WriteResultsCmgui.hpp"

//LatticeBoltzmann::LatticeBoltzmann()
//{
//  throw std::runtime_error("Declaring uninitialized LBM is not allowed");
//}

// cant use braces to initialize reference cuz gcc bug
// https://stackoverflow.com/questions/10509603/why-cant-i-initialize-a-
// reference-in-an-initializer-list-with-uniform-initializ
LatticeBoltzmann::LatticeBoltzmann(double t_total
  , const std::vector<std::vector<std::size_t>> &obstacles_position
  , bool is_ns
  , bool is_cd
  , bool is_taylor
  , bool is_instant
  , bool has_obstacles
  , LatticeModel &lm
  , CollisionNS &ns
  , CollisionCD &cd)
  : f {},
    boundary_f {},
    g {},
    boundary_g {},
    obstacles {},
    total_time_ {t_total},
    is_ns_ {is_ns},
    is_cd_ {is_cd},
    is_taylor_ {is_taylor},
    is_instant_ {is_instant},
    has_obstacles_ {has_obstacles},
    lm_ (lm),
    ns_ (ns),
    cd_ (cd)
{
  if (CheckParameters()) {
    throw std::runtime_error("Zero value in input parameters");
  }
  else {
    // Initializing variables with initial values
    if (is_ns_) f = ns_.lattice_eq;
    if (is_cd_) g = cd_.lattice_eq;
    if (has_obstacles) LatticeBoltzmann::Init(obstacles, obstacles_position);
  }
}

void LatticeBoltzmann::Init(std::vector<bool> &lattice
  , const std::vector<std::vector<std::size_t>> &position)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nd = lm_.GetNumberOfDimensions();
  obstacles.assign(nx * ny, false);
  if(position[0].size() != nd)
      throw std::runtime_error("Insufficient position information");
  for (auto pos : position) {
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    lattice[pos[1] * nx + pos[0]] = true;
  }  // pos
}

void LatticeBoltzmann::Obstacles(
    std::vector<std::vector<double>> &lattice)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (auto n = 0u; n < nx * ny; ++n) {
    if (obstacles[n]) {
      bool top = n / nx == ny - 1;
      bool bottom = n / nx == 0;
      bool left = n % nx == 0;
      bool right = n % nx == nx - 1;
      if (!top && !obstacles[n + nx]) lattice[n][N] = lattice[n + nx][S];
      if (!right && !obstacles[n + 1]) lattice[n][E] = lattice[n + 1][W];
      if (!left && !obstacles[n - 1]) lattice[n][W] = lattice[n - 1][E];
      if (!bottom && !obstacles[n - nx]) lattice[n][S] = lattice[n - nx][N];
      if (!top && !right && !obstacles[n + nx + 1])
          lattice[n][NE] = lattice[n + nx + 1][SW];
      if (!top && !left && !obstacles[n + nx - 1])
          lattice[n][NW] = lattice[n + nx - 1][SE];
      if (!bottom && !left && !obstacles[n - nx - 1])
          lattice[n][SW] = lattice[n - nx - 1][NE];
      if (!bottom && !right && !obstacles[n - nx + 1])
          lattice[n][SE] = lattice[n - nx + 1][NW];
    }
  }  // n
}

std::vector<std::vector<double>> LatticeBoltzmann::BoundaryCondition(
    const std::vector<std::vector<double>> &lattice)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nc = lm_.GetNumberOfDirections();
  std::vector<double> length_q(nc, 0.0);
  std::vector<std::vector<double>> left_boundary(ny, length_q);
  std::vector<std::vector<double>> top_boundary(nx, length_q);
  std::vector<std::vector<double>> corner_boundary(4, length_q);
  // initialize some boundaries to mirror their counterparts
  auto right_boundary(left_boundary);
  auto bottom_boundary(top_boundary);
  // no-slip boundary condition on left and right for taylor vortex
  // analytical solution, this is a temporary workaround
  if (is_taylor_) {
    for (auto y = 0u; y < ny; ++y) {
      auto n = y * nx;
      left_boundary[y][E] = lattice[n + 1][W];
      right_boundary[y][W] = lattice[n + nx - 1][E];
      if (y == 0) {
        left_boundary[y][SE] = lattice[nx * (ny - 1)][NW];
        right_boundary[y][SW] = lattice[nx * ny - 1][NE];
      }
      else {
        left_boundary[y][SE] = lattice[n - nx + 1][NW];
        right_boundary[y][SW] = lattice[n - 1][NE];
      }
      if (y == ny - 1) {
        left_boundary[y][NE] = lattice[0][SW];
        right_boundary[y][NW] = lattice[nx - 1][SE];
      }
      else {
        left_boundary[y][NE] = lattice[n + nx + 1][SW];
        right_boundary[y][NW] = lattice[n + nx + nx - 1][SE];
      }
    }  // y
  }
  // Periodic boundary condition on left and right
  else {
    for (auto y = 0u; y < ny; ++y) {
      auto n = y * nx;
      left_boundary[y][E] = lattice[n + nx - 1][E];
      left_boundary[y][NE] = lattice[n + nx - 1][NE];
      left_boundary[y][SE] = lattice[n + nx - 1][SE];
      right_boundary[y][W] = lattice[n][W];
      right_boundary[y][NW] = lattice[n][NW];
      right_boundary[y][SW] = lattice[n][SW];
    }  // y
  }
  // no-slip boundary condition on top and bottom
  for (auto x = 0u; x < nx; ++x) {
    auto n = (ny - 1) * nx;
    top_boundary[x][S] = lattice[x + n][N];
    bottom_boundary[x][N] = lattice[x][S];
    if (x == 0) {
      top_boundary[x][SW] = lattice[nx * ny - 1][NE];
      bottom_boundary[x][NW] = lattice[nx - 1][SE];
    }
    else {
      top_boundary[x][SW] = lattice[x + n - 1][NE];
      bottom_boundary[x][NW] = lattice[x - 1][SE];
    }
    if (x == nx - 1) {
      top_boundary[x][SE] = lattice[n][NW];
      bottom_boundary[x][NE] = lattice[0][SW];
    }
    else {
      top_boundary[x][SE] = lattice[x + n + 1][NW];
      bottom_boundary[x][NE] = lattice[x + 1][SW];
    }
  }  // x
  // mixture of boundaries at the corners
  corner_boundary[0][NE] = lattice[0][SW];
  corner_boundary[1][NW] = lattice[nx - 1][SE];
  corner_boundary[2][SE] = lattice[(ny - 1) * nx][NW];
  corner_boundary[3][SW] = lattice[ny * nx - 1][NE];

  // append the boundaries of different edges to the main boundary vector
  auto boundary(left_boundary);
  boundary.insert(end(boundary), begin(right_boundary), end(right_boundary));
  boundary.insert(end(boundary), begin(top_boundary), end(top_boundary));
  boundary.insert(end(boundary), begin(bottom_boundary), end(bottom_boundary));
  boundary.insert(end(boundary), begin(corner_boundary), end(corner_boundary));
  return boundary;
}

std::vector<std::vector<double>> LatticeBoltzmann::Stream(
    const std::vector<std::vector<double>> &lattice
  , const std::vector<std::vector<double>> &boundary)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nc = lm_.GetNumberOfDirections();
  // boundary indexes (left and right omitted)
  auto top = 2 * ny;
  auto bottom = 2 * ny + nx;
  auto corner = 2 * ny + 2 * nx;
  auto nx_big = nx + 2;  // number of columns in lattice with boundary
  auto ny_big = ny + 2;  // number of rows in lattice with boundary
  std::vector<std::vector<double>> lattice_with_bdr(nx_big * ny_big,
      std::vector<double>(nc, 0.0));
  // assemble lattice and boundary into a lattice with boundary
  // surrounding it
  for (auto y = 1u; y < ny_big - 1; ++y) {
    // left and right boundary
    lattice_with_bdr[y * nx_big] = boundary[y - 1];
    lattice_with_bdr[y * nx_big + nx + 1] = boundary[ny + y - 1];
    // main lattice
    for (auto x = 1u; x < nx_big - 1; ++x) {
      auto n = y * nx_big + x;
      auto m = (y - 1) * nx + x - 1;
      lattice_with_bdr[n] = lattice[m];
      // only fill in top and bottom boundary once, does it on the first
      // iteration of y loop
      if (y == 1) {
        // top and bottom boundary
        lattice_with_bdr[(ny_big - 1) * nx_big + x] = boundary[top + x - 1];
        lattice_with_bdr[x] = boundary[bottom + x - 1];
      }
    }  // x
  }  // y
  // corner boundary
  lattice_with_bdr[0] = boundary[corner];
  lattice_with_bdr[nx + 1] = boundary[corner + 1];
  lattice_with_bdr[(ny_big - 1) * nx_big] = boundary[corner + 2];
  lattice_with_bdr[ny_big * nx_big - 1] = boundary[corner + 3];
  // streams into temp lattice
  auto temp_lattice(lattice);
  for (auto n = 0u; n < nx * ny; ++n) {
    auto m = nx_big + 1 + n / nx * nx_big + n % nx;
    temp_lattice[n][E] = lattice_with_bdr[m - 1][E];
    temp_lattice[n][N] = lattice_with_bdr[m - nx_big][N];
    temp_lattice[n][W] = lattice_with_bdr[m + 1][W];
    temp_lattice[n][S] = lattice_with_bdr[m + nx_big][S];
    temp_lattice[n][NE] = lattice_with_bdr[m - nx_big - 1][NE];
    temp_lattice[n][NW] = lattice_with_bdr[m - nx_big + 1][NW];
    temp_lattice[n][SW] = lattice_with_bdr[m + nx_big + 1][SW];
    temp_lattice[n][SE] = lattice_with_bdr[m + nx_big - 1][SE];
  }  // n
  return temp_lattice;
}

void LatticeBoltzmann::TakeStep()
{
  if (is_ns_) {
    ns_.Collide(f);
    if (has_obstacles_) LatticeBoltzmann::Obstacles(f);
    f = LatticeBoltzmann::Stream(f, LatticeBoltzmann::BoundaryCondition(f));
    ns_.rho = lm_.ComputeRho(f);
    lm_.u = lm_.ComputeU(f, ns_.rho, ns_.source);
    ns_.ComputeEq();
  }
  if (is_cd_) {
    cd_.Collide(g);
    if (is_instant_) cd_.KillSource();
    if (has_obstacles_) LatticeBoltzmann::Obstacles(g);
    g = LatticeBoltzmann::Stream(g, LatticeBoltzmann::BoundaryCondition(g));
    cd_.rho = lm_.ComputeRho(g);
    cd_.ComputeEq();
  }
}

void LatticeBoltzmann::RunSim()
{
  auto dt = lm_.GetTimeStep();
  for (auto t = 0.0; t < total_time_; t += dt) LatticeBoltzmann::TakeStep();
}

void LatticeBoltzmann::RunSim(std::vector<std::vector<double>> &lattice)
{
  auto dt = lm_.GetTimeStep();
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto t_count = 0u;
  WriteResultsCmgui(lattice, nx, ny, t_count);
  for (auto t = 0.0; t < total_time_; t += dt) {
    LatticeBoltzmann::TakeStep();
    if (std::fmod(t, 0.002) < 1e-3) {
      WriteResultsCmgui(lattice, nx, ny, ++t_count);
      std::cout << t_count << " " << t << std::endl;
    }
  }
}

bool LatticeBoltzmann::CheckParameters()
{
  return total_time_ < 1e-20 || !(is_ns_ || is_cd_);
}
