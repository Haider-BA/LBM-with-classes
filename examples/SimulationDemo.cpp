#include <iostream>
#include "Lattice.hpp"
#include "UnitTest++.h"

SUITE(SimulationDemo)
{
static const std::size_t g_num_dimensions = 2;
static const std::size_t g_num_discrete_velocities = 9;
static const std::size_t g_num_rows = 21;
static const std::size_t g_num_cols = 31;
static const double g_dx = 0.0316;
static const double g_dt = 0.001;
static const double g_t_total = 1.0;
static const double g_diffusion_coefficient = 0.2;
static const double g_kinematic_viscosity = 0.2;
static const std::vector<double> g_u0 = {12, 0.2};
static const std::vector<double> g_u0_zero = {0, 0};
static const std::vector<std::vector<unsigned>> g_src_pos_f = {{0, 0}};
static const std::vector<std::vector<unsigned>> g_src_pos_g = {{15, 10}};
static const std::vector<std::vector<double>> g_src_strength_f = {{1, 1}};
static const std::vector<double> g_src_strength_g = {500};
static const bool g_is_cd = true;
static const bool g_is_ns = true;
static const bool g_is_instant = true;
static const bool g_is_not_cd = false;
static const bool g_is_not_ns = false;
static const bool g_is_not_instant = false;
static const double g_density_g = 1.0;
static const double g_density_f = 1.0;

TEST(DiffusionEquation)
{
  Lattice lattice(g_num_dimensions, g_num_discrete_velocities,
      g_num_rows, g_num_cols, g_dx, g_dt, g_t_total, g_diffusion_coefficient,
      g_kinematic_viscosity, g_density_f, g_density_g, g_u0_zero, g_src_pos_f,
      g_src_pos_g, g_src_strength_f, g_src_strength_g, g_is_cd, g_is_not_ns,
      g_is_not_instant);
  lattice.RunSim(lattice.g);
}

TEST(CoupledNSCDE)
{
  std::vector<std::vector<unsigned>> src_pos_f((g_num_cols * g_num_rows),
      {0, 0}) ;
  std::vector<std::vector<double>> src_strength_f((g_num_cols * g_num_rows),
      {10., -10.});
  for (auto y = 0u; y < g_num_rows; ++y) {
    for (auto x = 0u; x < g_num_cols; ++x) {
      auto n = y * g_num_cols + x;
      src_pos_f[n] = {x, y};
//      src_strength_f[n] = {(x % 5) * 5., (y % 5) * -5};
    }
  }
  Lattice lattice(g_num_dimensions, g_num_discrete_velocities,
      g_num_rows, g_num_cols, g_dx, g_dt, g_t_total, g_diffusion_coefficient,
      g_kinematic_viscosity, g_density_f, g_density_g, g_u0_zero, src_pos_f,
      g_src_pos_g, src_strength_f, g_src_strength_g, g_is_cd, g_is_ns,
      g_is_not_instant);
  lattice.RunSim(lattice.g);
}
TEST(NavierStokesEquation)
{
  std::vector<std::vector<unsigned>> src_pos_f((g_num_cols * g_num_rows),
      {0, 0}) ;
  std::vector<std::vector<double>> src_strength_f((g_num_cols * g_num_rows),
      {10., 0.});
  for (auto y = 0u; y < g_num_rows; ++y) {
    for (auto x = 0u; x < g_num_cols; ++x) {
      auto n = y * g_num_cols + x;
      src_pos_f[n] = {x, y};
//      src_strength_f[n] = {(x % 5) * 5., (y % 5) * -5};
    }
  }
  Lattice lattice(g_num_dimensions, g_num_discrete_velocities,
      g_num_rows, g_num_cols, g_dx, g_dt, g_t_total, g_diffusion_coefficient,
      g_kinematic_viscosity, g_density_f, g_density_g, g_u0_zero, src_pos_f,
      g_src_pos_g, src_strength_f, g_src_strength_g, g_is_not_cd, g_is_ns,
      g_is_not_instant);
  lattice.RunSim(lattice.u);
}

TEST(DiffusionEquationWithDiffusion)
{
  std::vector<std::vector<unsigned>> obs_pos;
  for (auto y = 0u; y < g_num_rows; ++y) {
    obs_pos.push_back({10, y});
    obs_pos.push_back({20, y});
  }
  for (auto x = 0u; x < g_num_cols; ++x) {
    obs_pos.push_back({x, 5});
  }
  double t_total = 0.002;
  Lattice lattice(g_num_dimensions, g_num_discrete_velocities,
      g_num_rows, g_num_cols, g_dx, g_dt, g_t_total, g_diffusion_coefficient,
      g_kinematic_viscosity, g_density_f, g_density_g, g_u0_zero, g_src_pos_f,
      g_src_pos_g, g_src_strength_f, g_src_strength_g, obs_pos, g_is_cd,
      g_is_not_ns, g_is_not_instant);

  lattice.RunSim(lattice.g);
//  lattice.RunSim();
}
}
