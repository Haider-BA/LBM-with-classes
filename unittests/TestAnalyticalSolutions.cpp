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
    , g_dt);
  CollisionNS ns(lm
    , src_pos_f
    , src_str_f
    , g_k_visco
    , g_rho0_f
    , u0);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g
    , u0);
  LatticeBoltzmann lbm(t_total
    , g_obs_pos
    , !g_is_ns
    , g_is_cd
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
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
}
