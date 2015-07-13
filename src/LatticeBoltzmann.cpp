#include "LatticeBoltzmann.hpp"
#include <cmath>  // std::fmod
#include <iomanip>  // std::setprecision
#include <iostream>
#include <stdexcept>  // std::runtime_error
#include <vector>
#include "Algorithm.hpp"
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
  , bool is_lid
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
    is_lid_ {is_lid},
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
  std::vector<std::vector<double>> left_bdr(ny, length_q);
  std::vector<std::vector<double>> top_bdr(nx, length_q);
  std::vector<std::vector<double>> corner_bdr(4, length_q);
  // initialize some boundaries to mirror their counterparts
  auto right_bdr(left_bdr);
  auto bottom_bdr(top_bdr);
  // Periodic boundary condition on left and right
  for (auto y = 0u; y < ny; ++y) {
    auto n = y * nx;
    left_bdr[y][E] = lattice[n + nx - 1][E];
    left_bdr[y][NE] = lattice[n + nx - 1][NE];
    left_bdr[y][SE] = lattice[n + nx - 1][SE];
    right_bdr[y][W] = lattice[n][W];
    right_bdr[y][NW] = lattice[n][NW];
    right_bdr[y][SW] = lattice[n][SW];
  }  // y
  // periodic boundary condition on top and bottom for taylor vortex
  // analytical solution, this is a temporary workaround
  if (is_taylor_) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = (ny - 1) * nx;
      top_bdr[x][S] = lattice[x][S];
      top_bdr[x][SW] = lattice[x][SW];
      top_bdr[x][SE] = lattice[x][SE];
      bottom_bdr[x][N] = lattice[n + x][N];
      bottom_bdr[x][NW] = lattice[n + x][NW];
      bottom_bdr[x][NE] = lattice[n + x][NE];
    }  // y
    // mixture of boundaries (periodic) at the corners
    corner_bdr[0][NE] = lattice[ny * nx - 1][NE];
    corner_bdr[1][NW] = lattice[(ny - 1) * nx][NW];
    corner_bdr[2][SE] = lattice[nx - 1][SE];
    corner_bdr[3][SW] = lattice[0][SW];
  }
  else {
    // no-slip boundary condition on top and bottom
    for (auto x = 0u; x < nx; ++x) {
      auto n = (ny - 1) * nx;
      top_bdr[x][S] = lattice[x + n][N];
      bottom_bdr[x][N] = lattice[x][S];
      top_bdr[x][SW] = lattice[(x == 0) ? nx * ny - 1 : x + n - 1][NE];
      bottom_bdr[x][NW] =lattice[(x == 0) ? nx - 1 : x - 1][SE];
      top_bdr[x][SE] = lattice[(x == nx - 1) ? n : x + n + 1][NW];
      bottom_bdr[x][NE] = lattice[(x == nx - 1) ? 0 : x + 1][SW];
    }  // x
    // mixture of boundaries (bounce-back) at the corners
    corner_bdr[0][NE] = lattice[0][SW];
    corner_bdr[1][NW] = lattice[nx - 1][SE];
    corner_bdr[2][SE] = lattice[(ny - 1) * nx][NW];
    corner_bdr[3][SW] = lattice[ny * nx - 1][NE];
  }
  // append the boundaries of different edges to the main boundary vector
  auto boundary(left_bdr);
  boundary.insert(end(boundary), begin(right_bdr), end(right_bdr));
  boundary.insert(end(boundary), begin(top_bdr), end(top_bdr));
  boundary.insert(end(boundary), begin(bottom_bdr), end(bottom_bdr));
  boundary.insert(end(boundary), begin(corner_bdr), end(corner_bdr));
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

void LatticeBoltzmann::BoundaryAndStream(
    std::vector<std::vector<double>> &lattice)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nc = lm_.GetNumberOfDirections();
  // Zou-He velocity BC for lid-driven flow, left, right and bottom boundary
  // from "http://lbmworkshop.com/wp-content/uploads/2011/08/Straight_
  // boundaries.pdf"
  // zero velocity
  auto u_x = 0.0;
  auto u_y = 0.0;
  auto u_lid = 0.01;
  for (auto y = 1u; y < ny - 1; ++y) {
    auto left = y * nx;
    auto right = left + nx - 1;
    auto eq_diff_left = 0.5 * (lattice[left][N] - lattice[left][S]);
    auto eq_diff_right = 0.5 * (lattice[right][N] - lattice[right][S]);
    auto rho_left = (lattice[left][0] + lattice[left][N] + lattice[left][S] +
        2.0 * (lattice[left][W] + lattice[left][NW] + lattice[left][SW])) /
        (1.0 - u_x);
    auto rho_right = (lattice[right][0] + lattice[right][N] +
        lattice[right][S] + 2.0 * (lattice[right][E] + lattice[right][NE] +
        lattice[right][SE])) / (1.0 - u_x);
    lattice[left][E] = lattice[left][W] + 2.0 / 3.0 * rho_left * u_x;
    lattice[left][NE] = lattice[left][SW] - eq_diff_left + rho_left *
        (u_x / 6.0 + 0.5 * u_y);
    lattice[left][SE] = lattice[left][NW] + eq_diff_left + rho_left *
        (u_x / 6.0 - 0.5 * u_y);
    lattice[right][W] = lattice[right][E] - 2.0 / 3.0 * rho_right * u_x;
    lattice[right][NW] = lattice[right][SE] - eq_diff_right - rho_right *
        (u_x / 6.0 - 0.5 * u_y);
    lattice[right][SW] = lattice[right][NE] + eq_diff_right - rho_right *
        (u_x / 6.0 + 0.5 * u_y);
  }
  for (auto x = 1u; x < nx - 1; ++x) {
    auto top = (ny - 1) * nx + x;
    auto eq_diff_bottom = 0.5 * (lattice[x][E] - lattice[x][W]);
    auto eq_diff_top = 0.5 * (lattice[top][E] - lattice[top][W]);
    auto rho_bottom = (lattice[x][0] + lattice[x][E] + lattice[x][W] + 2.0 *
        (lattice[x][S] + lattice[x][SW] + lattice[x][SE])) / (1 - u_y);
    auto rho_top = (lattice[top][0] + lattice[top][E] + lattice[top][W] + 2.0 *
        (lattice[top][S] + lattice[top][SW] + lattice[top][SE])) / (1 - u_y);
    lattice[x][N] = lattice[x][S] + 2.0 / 3.0 * rho_bottom * u_y;
    lattice[x][NE] = lattice[x][SW] - eq_diff_bottom + rho_bottom * (0.5 * u_x +
        u_y / 6.0);
    lattice[x][NW] = lattice[x][SE] + eq_diff_bottom - rho_bottom * (0.5 * u_x -
        u_y / 6.0);
    lattice[top][S] = lattice[top][N] - 2.0 / 3.0 * rho_top * u_y;
    lattice[top][SE] = lattice[top][NW] - eq_diff_top + rho_top * (0.5 * u_lid -
        u_y / 6.0);
    lattice[top][SW] = lattice[top][NE] + eq_diff_top - rho_top * (0.5 * u_lid +
        u_y / 6.0);
  }
  // Zou-He corners from "http://lbmworkshop.com/wp-content/uploads/2011/08/
  // Corners.pdf"
  // bottom-left corner
  lattice[0][E] = lattice[0][W] + 2.0 / 3.0 * u_x;
  lattice[0][N] = lattice[0][S] + 2.0 / 3.0 * u_y;
  lattice[0][NE] = lattice[0][SW] + (u_x + u_y) / 6.0;
  lattice[0][NW] = (u_y - u_x) / 12.0;
  lattice[0][SE] = (u_x - u_y) / 12.0;
  auto rho_1 = GetZerothMoment(lattice[nx + 1]);
  auto rho_2 = GetZerothMoment(lattice[2 * nx + 2]);
  // extrapolating density based on nodes along the diagonal
  auto rho_extrapolate = rho_1 - (rho_2 - rho_1);
  for (auto i = 1u; i < nc; ++i) rho_extrapolate -= lattice[0][i];
  lattice[0][0] = rho_extrapolate;
  // bottom-right corner
  lattice[nx - 1][W] = lattice[nx - 1][E] - 2.0 / 3.0 * u_x;
  lattice[nx - 1][N] = lattice[nx - 1][S] + 2.0 / 3.0 * u_y;
  lattice[nx - 1][NW] = lattice[nx - 1][SE] + (u_y - u_x) / 6.0;
  lattice[nx - 1][NE] = (u_x + u_y) / 12.0;
  lattice[nx - 1][SW] = (u_x + u_y) / -12.0;
  rho_1 = GetZerothMoment(lattice[2 * nx - 2]);
  rho_2 = GetZerothMoment(lattice[3 * nx - 3]);
  // extrapolating density based on nodes along the diagonal
  rho_extrapolate = rho_1 - (rho_2 - rho_1);
  for (auto i = 1u; i < nc; ++i) rho_extrapolate -= lattice[nx - 1][i];
  lattice[nx - 1][0] = rho_extrapolate;
  // top-left corner
  auto top_left = (ny - 1) * nx;
  lattice[top_left][E] = lattice[top_left][W] + 2.0 / 3.0 * u_lid;
  lattice[top_left][S] = lattice[top_left][N] - 2.0 / 3.0 * u_y;
  lattice[top_left][SE] = lattice[top_left][NW] + (u_lid - u_y) / 6.0;
  lattice[top_left][NE] = (u_lid + u_y) / 12.0;
  lattice[top_left][SW] = (u_lid + u_y) / -12.0;
  rho_1 = GetZerothMoment(lattice[top_left - nx + 1]);
  rho_2 = GetZerothMoment(lattice[top_left - 2 * nx + 2]);
  // extrapolating density based on nodes along the diagonal
  rho_extrapolate = rho_1 - (rho_2 - rho_1);
  for (auto i = 1u; i < nc; ++i) rho_extrapolate -= lattice[top_left][i];
  lattice[top_left][0] = rho_extrapolate;
  // top-right corner
  auto top_right = nx * ny - 1;
  lattice[top_right][W] = lattice[top_right][E] - 2.0 / 3.0 * u_lid;
  lattice[top_right][S] = lattice[top_right][N] - 2.0 / 3.0 * u_y;
  lattice[top_right][SW] = lattice[top_right][NE] - (u_lid + u_y) / 6.0;
  lattice[top_right][NW] = (u_y - u_lid) / 12.0;
  lattice[top_right][SE] = (u_lid - u_y) / 12.0;
  rho_1 = GetZerothMoment(lattice[top_right - nx - 1]);
  rho_2 = GetZerothMoment(lattice[top_right - 2 * nx - 2]);
  // extrapolating density based on nodes along the diagonal
  rho_extrapolate = rho_1 - (rho_2 - rho_1);
  for (auto i = 1u; i < nc; ++i) rho_extrapolate -= lattice[top_right][i];
  lattice[top_right][0] = rho_extrapolate;
  // Streaming
}

void LatticeBoltzmann::TakeStep()
{
  if (is_ns_) {
    ns_.Collide(f);
    if (has_obstacles_) LatticeBoltzmann::Obstacles(f);
    if (is_lid_) {
      LatticeBoltzmann::BoundaryAndStream(f);
    }
    else {
      f = LatticeBoltzmann::Stream(f, LatticeBoltzmann::BoundaryCondition(f));
    }
    ns_.rho = lm_.ComputeRho(f);
    lm_.u = lm_.ComputeU(f, ns_.rho, ns_.source);
    ns_.ComputeEq();
  }
  if (is_cd_) {
    cd_.Collide(g);
    if (is_instant_) cd_.KillSource();
    if (has_obstacles_) LatticeBoltzmann::Obstacles(g);
    if (is_lid_) {
      LatticeBoltzmann::BoundaryAndStream(g);
    }
    else {
      g = LatticeBoltzmann::Stream(g, LatticeBoltzmann::BoundaryCondition(g));
    }
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
