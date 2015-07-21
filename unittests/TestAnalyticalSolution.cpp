#include <cmath>  // exp, sin, cos
#include <iomanip>
#include <iostream>
#include <vector>
#include "BounceBackNodes.hpp"
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "CollisionNSF.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "StreamPeriodic.hpp"
#include "UnitTest++.h"
#include "WriteResultsCmgui.hpp"

SUITE(TestAnalyticalSolutions)
{
static const auto g_dx = 0.0316;
static const auto g_dt = 0.001;
static const auto g_cs_sqr = (g_dx / g_dt) * (g_dx / g_dt) / 3.0;
static const auto g_k_visco = 0.2;
static const auto g_rho0_f = 1.0;
static const auto g_d_coeff = 0.2;
static const auto g_rho0_g = 1.0;
static const auto g_is_instant = true;
static const auto g_is_prestream = true;
static const auto g_pi = 3.14159265;
static const auto g_2pi = 3.14159265 * 2;

TEST(AnalyticalDiffusion)
{
  std::size_t ny = 201;
  std::size_t nx = 201;
  std::vector<std::vector<std::size_t>> src_pos_g = {{100, 100}};
  std::vector<double> src_str_g = {1000.0};  // unit conversion
  std::vector<double> u0 = {0.0, 0.0};
  auto time_steps = 100;
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamPeriodic sp(lm);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g
    , g_is_instant);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  // ignoring comparison with t = 0 since it will cause divide by zero error
  g.TakeStep();
  for (auto t = 1; t < time_steps; ++t) {
    g.TakeStep();
    auto n = 0;
    std::cout << t << std::endl;
    for (auto node : cd.rho) {
      auto y = abs(n / nx - 100);
      auto x = abs(n % nx - 100);
      // Analytical solution from http://nptel.ac.in/courses/105103026/34
      double rho = exp(-1.0 * (y * y + x * x) / 4.0 / g_d_coeff / t) /
          (4.0 * g_pi * t * g_d_coeff);
      // ignoring the first 6 time steps
      if (t > 6) CHECK_CLOSE(rho, node - 1.0, 0.0068);
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
  std::vector<double> u0 = {0, 0};
  auto time_steps = 1000;
  for (auto n = 0u; n < nx * ny; ++n) src_pos_f.push_back({n % nx, n / nx});
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamPeriodic sp(lm);
  CollisionNSF nsf(lm
    , src_pos_f
    , src_str_f
    , g_k_visco
    , g_rho0_f);
  BounceBackNodes bbnsf(g_is_prestream
    , nsf
    , lm);
  LatticeBoltzmann f(lm
    , nsf
    , sp);
  for (auto x = 0u; x < nx; ++x) {
    bbnsf.AddNode(x, 0);
    bbnsf.AddNode(x, ny - 1);
  }
  f.AddBoundaryNodes(&bbnsf);
  for (auto t = 0; t < 2000; ++t) {
    f.TakeStep();
  }
  // calculation of analytical u_max according to formula in Guo2002 after
  auto length = static_cast<double>(ny / 2) - 1.0;
  double u_max = body_force / 1000 * length * length / 2 / g_k_visco;
  for (auto x = 0u; x < 1; ++x) {
    double y_an = 0.5;
    for (auto y = 1u; y < ny - 1; ++y) {
      auto n = y * nx + x;
      double u_an = u_max * (1.0 - (y_an - length) * (y_an - length) /
          (length * length));
      auto u_sim = lm.u[n][0] + lm.u[n][1];
      CHECK_CLOSE(u_an, u_sim, u_an * 0.03);
      if (y < 8) ++y_an;
      if (y > 8) --y_an;
    }  // y
  }  // x
}
}
