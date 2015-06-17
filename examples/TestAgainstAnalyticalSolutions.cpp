#include <cmath>  // exp
#include <iostream>
#include "Lattice.hpp"
#include "UnitTest++.h"

SUITE(AnalyticalSolutionTests)
{
static const double zero_tol = 1e-20;
static const double loose_tol = 1e-5;
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
static const std::size_t g_num_dimensions = 2;
static const std::size_t g_num_discrete_velocities = 9;
static const std::size_t g_num_rows = 18  ;
static const std::size_t g_num_cols = 34;
static const double g_dx = 0.0316;
static const double g_dt = 0.001;
static const double g_t_total = 1.0;
static const double g_diffusion_coefficient = 0.2;
static const double g_kinematic_viscosity = 0.2;
static const std::vector<double> g_u0 = {123, 321};
static const std::vector<double> g_u0_zero = {0, 0};
static const std::vector<std::vector<unsigned>> g_src_pos_f = {{0, 0}};
static const std::vector<std::vector<unsigned>> g_src_pos_g = {{0, 0}};
static const std::vector<std::vector<double>> g_src_strength_f = {{1, 1}};
static const std::vector<double> g_src_strength_g = {1000};
static const bool g_is_cd = true;
static const bool g_is_ns = true;
static const bool g_is_instant = true;
static const bool g_is_not_cd = false;
static const bool g_is_not_ns = false;
static const bool g_is_not_instant = false;
static const double g_density_g = 2.0;
static const double g_density_f = 2.0;
static const Lattice g_lattice(g_num_dimensions, g_num_discrete_velocities,
      g_num_rows, g_num_cols, g_dx, g_dt, g_t_total, g_diffusion_coefficient,
      g_kinematic_viscosity, g_density_f, g_density_g, g_u0, g_src_pos_f,
      g_src_pos_g, g_src_strength_f, g_src_strength_g, g_is_cd, g_is_ns,
      g_is_not_instant);

TEST(DiffusionEquation)
{
  std::size_t nx = 201;
  std::size_t ny = 201;
  double density_g = 1.0;
  std::vector<std::vector<unsigned>> src_pos_g = {{100, 100}};
  Lattice lattice(g_num_dimensions, g_num_discrete_velocities, ny, nx, g_dx,
      g_dt, g_t_total, g_diffusion_coefficient, g_kinematic_viscosity,
      g_density_f, density_g, g_u0_zero, g_src_pos_f, src_pos_g,
      g_src_strength_f, g_src_strength_g, g_is_cd, g_is_not_ns, g_is_instant);
  lattice.InitAll();
  for (auto t = 0; t < 100; ++t) {
    lattice.TakeStep();
    if (t) {
      for (auto y = 0u; y < ny; ++y) {
        for (auto x = 0u; x < nx; ++x) {
          auto n = y * nx + x;
          auto y_an = abs(y - 100);
          auto x_an = abs(x - 100);
          // Analytical solution from http://nptel.ac.in/courses/105103026/34
          double rho_an = 1.0 + exp(-1.0 * (y_an * y_an + x_an * x_an) / 4.0 /
               g_diffusion_coefficient / t) / (4.0 * 3.1415926 * t *
               g_diffusion_coefficient);
          if(!(y == 100 && x == 100) && t > 3) {
            CHECK_CLOSE(rho_an, lattice.rho_g[n], 0.007);
          }
        }  // x
      }  // y
    }
  }  // t
}

TEST(PoiseiuilleFlow)
{
  double body_force = 1.0;
  std::vector<std::vector<unsigned>> src_pos_f((g_num_cols * g_num_rows),
      {0, 0});
  std::vector<std::vector<double>> src_strength_f((g_num_cols * g_num_rows),
      {body_force, 0});
  for (auto y = 0u; y < g_num_rows; ++y) {
    for (auto x = 0u; x < g_num_cols; ++x) {
      auto n = y * g_num_cols + x;
      src_pos_f[n] = {x, y};
    }
  }
  double density_f = 1.0;
  Lattice lattice(g_num_dimensions, g_num_discrete_velocities, g_num_rows,
      g_num_cols, g_dx, g_dt, g_t_total, g_diffusion_coefficient,
      g_kinematic_viscosity, density_f, g_density_g, g_u0_zero, src_pos_f,
      g_src_pos_g, src_strength_f, g_src_strength_g, g_is_not_cd, g_is_ns,
      g_is_instant);
  lattice.RunSim();
  // calculation of analytical u_max?
  double u_max = body_force / 1000 * (g_num_rows / 2) * (g_num_rows / 2) / 2 /
      g_kinematic_viscosity;
  for (auto x = 0u; x < g_num_cols; ++x) {
    unsigned y_an = 1;
    for (auto y = 0u; y < g_num_rows; ++y) {
      auto n = y * g_num_cols + x;
      double u_an = u_max * (1.0 - (y_an - 10.0) * (y_an - 10.0) / 100.0);
      double u_sim = lattice.u[n][0] + lattice.u[n][1];
      // accuracy of answer?
      CHECK_CLOSE(u_an, u_sim, 0.02);
      if (y < 8) ++y_an;
      if (y > 8) --y_an;
    }  // y
  }  // x
}
}
