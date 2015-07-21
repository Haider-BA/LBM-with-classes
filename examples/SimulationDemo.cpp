#include <iostream>
#include <vector>
#include "BounceBackNodes.hpp"
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "CollisionNSF.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "Printing.hpp"
#include "StreamPeriodic.hpp"
#include "UnitTest++.h"
#include "WriteResultsCmgui.hpp"

SUITE(SimulationDemo)
{
const static auto g_dx = 0.0316;
const static auto g_dt = 0.001;
const static auto g_k_visco = 0.2;
const static auto g_rho0_f = 1.0;
const static auto g_d_coeff = 0.2;
const static auto g_rho0_g = 1.0;
const static auto g_time_steps = 1000;
const static auto g_is_instant = true;
const static auto g_is_prestream = true;

TEST(SimulateDiffusion)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_g = {{15, 10}};
  std::vector<double> src_str_g = {50};
  std::vector<double> u0 = {0, 0};
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
    , !g_is_instant);
  BounceBackNodes bbcd(g_is_prestream
    , cd
    , lm);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  for (auto t = 0u; t < 501; ++t) {
    g.TakeStep();
    WriteResultsCmgui(g.df, nx, ny, t);
  }
}

TEST(SimulateConvectionDiffusion)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_g = {{15, 10}};
  std::vector<double> src_str_g = {50};
  std::vector<double> u0 = {10, 10};
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
    , !g_is_instant);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  for (auto t = 0u; t < 501; ++t) {
    g.TakeStep();
    WriteResultsCmgui(g.df, nx, ny, t);
  }
}

TEST(SimulatePoiseuilleFlow)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {10.0, 0.0});
  std::vector<double> u0 = {0.0, 0.0};
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
  for (auto t = 0u; t < 501; ++t) {
    f.TakeStep();
    WriteResultsCmgui(lm.u, nx, ny, t);
  }
}

TEST(SimulateNSCDCoupling)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {50.0, -10.0});
  std::vector<std::vector<std::size_t>> src_pos_g = {{15, 10}};
  std::vector<double> src_str_g = {50};
  std::vector<double> u0 = {-5.0, 0.0};
  for (auto n = 0u; n < nx * ny; ++n) {
    src_pos_f.push_back({n % nx, static_cast<std::size_t>(n / nx)});
  }
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
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g
    , !g_is_instant);
  BounceBackNodes bbnsf(g_is_prestream
    , nsf
    , lm);
  BounceBackNodes bbcd(g_is_prestream
    , cd
    , lm);
  LatticeBoltzmann f(lm
    , nsf
    , sp);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  for (auto x = 0u; x < nx; ++x) {
    bbnsf.AddNode(x, 0);
    bbnsf.AddNode(x, ny - 1);
    bbcd.AddNode(x, 0);
    bbcd.AddNode(x, ny - 1);
  }
  f.AddBoundaryNodes(&bbnsf);
  g.AddBoundaryNodes(&bbcd);
  for (auto t = 0u; t < 501; ++t) {
    f.TakeStep();
    g.TakeStep();
    WriteResultsCmgui(g.df, nx, ny, t);
  }
}
}
