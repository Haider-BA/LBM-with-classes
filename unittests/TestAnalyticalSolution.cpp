#include <cmath>  // exp, sin, cos
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "BouncebackNodes.hpp"
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "CollisionNSF.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "StreamD2Q9.hpp"
#include "StreamPeriodic.hpp"
#include "UnitTest++.h"
#include "WriteResultsCmgui.hpp"
#include "ZouHeNodes.hpp"

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
static const auto g_pi = 3.14159265358979323846;
static const auto g_2pi = g_pi * 2;

TEST(AnalyticalDiffusion)
{
  // checking at end of simulation, sum(abs(rho_sim - rho_an)) / sum(rho_an)
  // appears to be unable to display rho values, rho < 5e-6
  // dx value cannot be too small, error greatest at src
  std::size_t ny = 201;
  std::size_t nx = 201;
  auto dx = 0.0316;
  auto dt = dx * dx;
  auto d_coeff = 0.05;
  auto src_g_an = 1.0;
  auto src_g = src_g_an / dt / dx / dx;
  auto src_coord = static_cast<std::size_t>(nx / 2);
  std::vector<std::vector<std::size_t>> src_pos_g = {{src_coord, src_coord}};
  std::vector<double> src_str_g = {src_g};  // unit conversion
  std::vector<double> u0 = {0.0, 0.0};
  auto time_steps = 3000;
  LatticeD2Q9 lm(ny
    , nx
    , dx
    , dt
    , u0);
  StreamPeriodic sp(lm);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , d_coeff
    , g_rho0_g
    , g_is_instant);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  for (auto t = 0; t < time_steps; ++t) g.TakeStep();
  std::ofstream myfile;
  myfile.open("diffusion.csv");
  myfile << "rho_y,rho_x" << std::endl;
  for (auto i = 0u; i < nx; ++i) {
    auto y = src_coord * nx + i;
    auto x = i * nx + src_coord;
    myfile << cd.rho[y] << "," << cd.rho[x] << std::endl;
  }
  myfile.close();
  auto n = 0;
  auto std_error = 0.0;
  auto ana_sum = 0.0;
  auto t_an = static_cast<double>(time_steps) * g_dt;
  for (auto node : cd.rho) {
    auto y = abs(n / nx - src_coord);
    auto x = abs(n % nx - src_coord);
    auto y_an = static_cast<double>(y) * g_dx;
    auto x_an = static_cast<double>(x) * g_dx;
    // Analytical solution from http://nptel.ac.in/courses/105103026/34
    // since src_g is g/m^2/s
    double rho_an = src_g_an * exp(-1.0 * (y_an * y_an + x_an * x_an) / 4.0 /
        d_coeff / t_an) / (4.0 * g_pi * t_an * d_coeff);
    std_error += fabs(node - rho_an - 1.0);
    ana_sum += rho_an;
    ++n;
  }  // n
  CHECK_CLOSE(0.0, std_error / ana_sum, 1e-3);
}

TEST(AnalyticalPoiseuille)
{
  std::size_t ny = 18;
  std::size_t nx = 34;
  double body_force = 10.0;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {body_force, 0});
  std::vector<double> u0 = {0, 0};
  auto time_steps = 3000;
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
  BouncebackNodes bbnsf(lm
    , &sp);
  LatticeBoltzmann f(lm
    , nsf
    , sp);
  for (auto x = 0u; x < nx; ++x) {
    bbnsf.AddNode(x, 0);
    bbnsf.AddNode(x, ny - 1);
  }
  f.AddBoundaryNodes(&bbnsf);
  for (auto t = 0; t < time_steps; ++t) {
    f.TakeStep();
    if (t % 6 == 0) WriteResultsCmgui(lm.u, nx, ny, t / 6);
  }
  // calculation of analytical u_max according to formula in Guo2002
  auto length = static_cast<double>(ny / 2);
  auto length_an = static_cast<double>(ny / 2) * g_dx;
  auto visco_an = g_k_visco * g_dx * g_dx / g_dt;
  double u_max = body_force * length_an * length_an / 2 / visco_an;
  // check against velocities in the middle of the channel
  for (auto x = 10u; x < nx - 10; ++x) {
    for (auto y = 0u; y < ny; ++y) {
      auto n = y * nx + x;
      auto y_an = fabs(static_cast<double>(y) - length + 0.5) * g_dx;
      double u_an = u_max * (1.0 - y_an * y_an / (length_an * length_an));
      auto u_sim = lm.u[n][0] + lm.u[n][1];
      std::cout << u_sim << std::endl;
      CHECK_CLOSE(u_an, u_sim, u_an * 0.02);
    }  // y
  }  // x
}

TEST(AnalyticalPoiseuilleZH)
{
  std::size_t ny = 38;
  std::size_t nx = 150;
  std::vector<double> u0 = {0.0, 0.0};
  auto u_in = 1.35;
  auto time_steps = 3000;
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamD2Q9 sd(lm);
  StreamPeriodic sp(lm);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  BouncebackNodes fwbb(lm
    , &ns);
  ZouHeNodes inlet(lm
    , ns);
  ZouHeNodes outlet(lm
    , ns);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  for (auto x = 0u; x < nx; ++x) {
    fwbb.AddNode(x, 0);
    fwbb.AddNode(x, ny - 1);
  }
  for (auto y = 1u; y < ny - 1; ++y) {
    inlet.AddNode(0, y, u_in, 0);
    outlet.AddNode(nx - 1, y, 0.0, 0.0);
  }
  outlet.ToggleNormalFlow();
  f.AddBoundaryNodes(&fwbb);
  f.AddBoundaryNodes(&inlet);
  f.AddBoundaryNodes(&outlet);
  for (auto t = 0; t < time_steps; ++t) f.TakeStep();
  auto length = static_cast<double>(ny / 2 - 1);
  auto length_an = static_cast<double>(ny / 2 - 1) * g_dx;
  double u_max = u_in * 1.5;
  // check against velocities in the middle of the channel
  for (auto x = 100; x < 101; ++x) {
    for (auto y = 1u; y < ny - 1; ++y) {
      auto n = y * nx + x;
      auto y_an = fabs(static_cast<double>(y - 1) - length + 0.5) * g_dx;
      double u_an = u_max * (1.0 - y_an * y_an / (length_an * length_an));
      auto u_sim = lm.u[n][0];
      CHECK_CLOSE(u_an, u_sim, u_an * 0.025);
      std::cout << u_sim << std::endl;
    }  // y
  }  // x
}

TEST(AnalyticalTaylorVortex)
{
  // have to use odd number for sizes
  std::size_t ny = 65;
  std::size_t nx = 65;
  auto time_steps = 300;
  // analytical solution parameters
  std::vector<std::vector<double>> u_lattice_an;
  std::vector<double> rho_lattice_an;
  auto u0_an = 0.001;
  auto k_visco = 0.25;
  // using one k since it's a square box
  auto k = g_2pi / nx;
  for (auto n = 0u; n < nx * ny; ++n) {
    auto x_an = static_cast<double>(n % nx) * k;
    auto y_an = static_cast<double>(n / nx) * k;
    // analytical formula from "Interpolation methods and the accuracy of
    // lattice-Boltzmann mesh refinement" eq17
    auto u_an = -1.0 * u0_an * cos(x_an) * sin(y_an);
    auto v_an = u0_an * sin(x_an) * cos(y_an);
    u_lattice_an.push_back({u_an, v_an});
    auto rho_an = g_rho0_f - 0.25 / g_cs_sqr * u0_an * u0_an *
        (cos(2.0 * x_an) + cos(2.0 * y_an));
    rho_lattice_an.push_back(rho_an);
  }  // n
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u_lattice_an);
  StreamPeriodic sp(lm);
  CollisionNS ns(lm
    , k_visco
    , rho_lattice_an);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  // According to t_c formula in Guo2002 pg5, checks simulation results against
  // analytical result until velocity is 25% of initial value
  for (auto t = 0; t < time_steps; ++t) {
    f.TakeStep();
    if (t == 148 || t == 299) {
      for (auto y = 0u; y < ny; ++y) {
        auto n = y * nx;
        std::cout << lm.u[n][0] << std::endl;
      }
    }
    WriteResultsCmgui(lm.u, nx, ny, t);
    for (auto n = 0u; n < nx * ny; ++n) {
        auto x_an = static_cast<double>(n % nx) * k;
        auto y_an = static_cast<double>(n / nx) * k;
        auto u_an = -1.0 * u0_an * cos(x_an) * sin(y_an) *
            exp(-2.0 * k_visco * k * k * t);
        auto v_an = u0_an * sin(x_an) * cos(y_an) *
            exp(-2.0 * k_visco * k * k * t);
      // checks that simulation is within 1% of analytical value if analytical
      // value is not zero, else check they are less than 1e-8
      CHECK_CLOSE(u_an, lm.u[n][0], (fabs(u_an) > 1e-20) ? fabs(u_an) * 0.01 :
          1e-8);
      CHECK_CLOSE(v_an, lm.u[n][1], (fabs(v_an) > 1e-20) ? fabs(v_an) * 0.01 :
          1e-8);
    }  // y
  }  // t
}

TEST(AnalyticalTaylorVortexForce)
{
  // have to use odd number for sizes
  std::size_t ny = 65;
  std::size_t nx = 65;
  auto time_steps = 300;
  // analytical solution parameters
  std::vector<std::vector<double>> u_lattice_an;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f;
  auto u0_an = 0.001;
  auto body_force = u0_an * u0_an;
  auto k_visco = 0.25;
  auto rho0_f = 1000000.0;
  // using one k since it's a square box
  auto k = g_2pi / nx;
  for (auto n = 0u; n < nx * ny; ++n) {
    auto x = n % nx;
    auto y = n / nx;
    src_pos_f.push_back({x, y});
    auto x_an = static_cast<double>(x) * k;
    auto y_an = static_cast<double>(y) * k;
    // analytical formula from "Interpolation methods and the accuracy of
    // lattice-Boltzmann mesh refinement" eq17
    auto u_an = -1.0 * u0_an * cos(x_an) * sin(y_an);
    auto v_an = u0_an * sin(x_an) * cos(y_an);
    u_lattice_an.push_back({u_an, v_an});
    auto f_x = -0.5 / g_cs_sqr * body_force * sin(2.0 * x_an);
    auto f_y = -0.5 / g_cs_sqr * body_force * sin(2.0 * y_an);
    src_str_f.push_back({f_x, f_y});
  }  // n
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u_lattice_an);
  StreamPeriodic sp(lm);
  CollisionNSF nsf(lm
    , src_pos_f
    , src_str_f
    , k_visco
    , rho0_f);
  LatticeBoltzmann f(lm
    , nsf
    , sp);
  // According to t_c formula in Guo2002 pg5, checks simulation results against
  // analytical result until velocity is 25% of initial value
  for (auto t = 0; t < time_steps; ++t) {
    for (auto n = 0u; n < nx * ny; ++n) {
      auto x_an = static_cast<double>(n % nx) * k;
      auto y_an = static_cast<double>(n / nx) * k;
      // analytical formula from "Interpolation methods and the accuracy of
      // lattice-Boltzmann mesh refinement" eq17
      auto f_x = -0.5 / g_cs_sqr * body_force * sin(2.0 * x_an) *
          exp(-2.0 * k * k * k_visco * t);
      auto f_y = -0.5 / g_cs_sqr * body_force * sin(2.0 * y_an) *
          exp(-2.0 * k * k * k_visco * t);
      src_str_f[n] = {f_x, f_y};
    }  // n
    for (auto n = 0u; n < nx * ny; ++n) {
      nsf.source[n][0] /= nsf.rho[n];
      nsf.source[n][1] /= nsf.rho[n];
    }
    nsf.InitSource(src_pos_f
      , src_str_f);
    f.TakeStep();
    for (auto n = 0u; n < nx * ny; ++n) {
        auto x_an = static_cast<double>(n % nx) * k;
        auto y_an = static_cast<double>(n / nx) * k;
        auto u_an = -1.0 * u0_an * cos(x_an) * sin(y_an) *
            exp(-2.0 * k_visco * k * k * t);
        auto v_an = u0_an * sin(x_an) * cos(y_an) *
            exp(-2.0 * k_visco * k * k * t);
      // checks that simulation is within 1% of analytical value if analytical
      // value is not zero, else check they are less than 1e-7
      CHECK_CLOSE(u_an, lm.u[n][0], (fabs(u_an) > 1e-20) ? fabs(u_an) * 0.01 :
          1e-7);
      CHECK_CLOSE(v_an, lm.u[n][1], (fabs(v_an) > 1e-20) ? fabs(v_an) * 0.01 :
          1e-7);
    }  // y
  }  // t
}
}
