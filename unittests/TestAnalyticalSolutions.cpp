#include <cmath>  // exp
#include <iomanip>
#include <iostream>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"
#include "UnitTest++.h"

SUITE(TestAnalyticalSolutions)
{
const static double g_dx = 0.0316;
const static double g_dt = 0.001;
const static double g_k_visco = 0.2;
const static double g_rho0_f = 1.0;
const static double g_d_coeff = 0.2;
const static double g_rho0_g = 1.0;
const static std::vector<std::vector<std::size_t>> g_obs_pos;
static const bool g_is_ns = true;
static const bool g_is_cd = true;
static const bool g_is_instant = true;
static const bool g_no_obstacles = false;

TEST(AnalyticalDiffusion)
{
  std::size_t ny = 201;
  std::size_t nx = 201;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f;
  std::vector<std::vector<std::size_t>> src_pos_g{{100, 100}};
  std::vector<double> src_str_g{1000};  // unit conversion
  std::vector<double> u0{0, 0};
  double t_total{1.0};
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  CollisionNS ns(lm
    , src_pos_f
    , src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(t_total
    , g_obs_pos
    , !g_is_ns
    , g_is_cd
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  // ignoring comparison with t = 0 since it will cause divide by zero error
  lbm.TakeStep();
  for (auto t = 1; t < 100; ++t) {
    lbm.TakeStep();
    auto n = 0;
    for (auto node : cd.rho) {
      auto y = abs(n / nx - 100);
      auto x = abs(n % nx - 100);
      // Analytical solution from http://nptel.ac.in/courses/105103026/34
      double rho = exp(-1.0 * (y * y + x * x) / 4.0 / g_d_coeff / t) /
          (4.0 * 3.1415926 * t * g_d_coeff);
      if (t > 3) CHECK_CLOSE(rho, node - 1.0, 0.0068);
      ++n;
    }  // n
  }  // t
}

TEST(AnalyticalPoiseuille)
{
  std::size_t ny = 18;
  std::size_t nx = 34;
  double body_force = 10.0;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {body_force, 0});
  std::vector<std::vector<std::size_t>> src_pos_g;
  std::vector<double> src_str_g;
  std::vector<double> u0{0, 0};
  double t_total{1.0};
  for (auto n = 0u; n < nx * ny; ++n) src_pos_f.push_back({n % nx, n / nx});
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  CollisionNS ns(lm
    , src_pos_f
    , src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(t_total
    , g_obs_pos
    , g_is_ns
    , !g_is_cd
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.RunSim();
  // calculation of analytical u_max according to formula in Guo2002 after
  double u_max = body_force / 1000 * (ny / 2) * (ny / 2) / 2 / g_k_visco;
  for (auto x = 0u; x < nx; ++x) {
    double y_an = 0.5;
    for (auto y = 0u; y < ny; ++y) {
      auto n = y * nx + x;
      double u_an = u_max * (1.0 - (y_an - 9.0) * (y_an - 9.0) / 81.0);
      auto u_sim = lm.u[n][0] + lm.u[n][1];
      // simulation value is within 2.35% of analytical value
      CHECK_CLOSE(u_an, u_sim, u_an * 0.0235);
      if (y < 8) ++y_an;
      if (y > 8) --y_an;
    }  // y
  }  // x
}
}
