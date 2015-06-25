#include <iostream>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"
#include "UnitTest++.h"

SUITE(SimulationDemo)
{
const static double g_dx = 0.0316;
const static double g_dt = 0.001;
const static std::vector<std::vector<std::size_t>> g_src_pos_f;
const static std::vector<std::vector<double>> g_src_str_f;
const static double g_k_visco = 0.2;
const static double g_rho0_f = 1.0;
const static std::vector<std::vector<std::size_t>> g_src_pos_g;
const static std::vector<double> g_src_str_g;
const static double g_d_coeff = 0.2;
const static double g_rho0_g = 1.0;
const static double g_t_total = 1.0;
const static std::vector<std::vector<std::size_t>> g_obs_pos;
static const bool g_is_ns = true;
static const bool g_is_cd = true;
static const bool g_is_instant = true;
static const bool g_no_obstacles = false;

TEST(SimulateDiffusion)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_g{{15, 10}};
  std::vector<double> src_str_g{50};
  std::vector<double> u0{0, 0};
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_obs_pos
    , !g_is_ns
    , g_is_cd
    , !g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.RunSim(lbm.g);
}

TEST(SimulateConvectionDiffusion)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_g{{15, 10}};
  std::vector<double> src_str_g{50};
  std::vector<double> u0{10, -0.01};
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_obs_pos
    , !g_is_ns
    , g_is_cd
    , !g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.RunSim(lbm.g);
}

TEST(SimulatePoiseuilleFlow)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {10.0, 0.0});
  std::vector<double> u0{0.0, 0.0};
  for (auto n = 0u; n < nx * ny; ++n) {
    src_pos_f.push_back({n % nx, static_cast<std::size_t>(n / nx)});
  }
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
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_obs_pos
    , g_is_ns
    , !g_is_cd
    , !g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.RunSim(lm.u);
}

TEST(SimulateNSCDCoupling)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {50.0, -10.0});
  std::vector<std::vector<std::size_t>> src_pos_g{{15, 10}};
  std::vector<double> src_str_g{50};
  std::vector<double> u0{-5.0, 0.0};
  for (auto n = 0u; n < nx * ny; ++n) {
    src_pos_f.push_back({n % nx, static_cast<std::size_t>(n / nx)});
  }
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
  LatticeBoltzmann lbm(g_t_total
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.RunSim(lbm.g);
}

TEST(SimulateNSCDCouplingWithObstacles)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {50.0, -10.0});
  std::vector<std::vector<std::size_t>> src_pos_g{{15, 10}};
  std::vector<double> src_str_g{50};
  std::vector<std::vector<std::size_t>> obs_pos;
  std::vector<double> u0{-5.0, 0.0};
  for (auto n = 0u; n < nx * ny; ++n) {
    src_pos_f.push_back({n % nx, static_cast<std::size_t>(n / nx)});
  }
  for (auto y = 9u; y < 12u; ++y) {
    obs_pos.push_back({17, y});
    obs_pos.push_back({18, y});
    obs_pos.push_back({19, y});
  }
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
  LatticeBoltzmann lbm(g_t_total
    , obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_instant
    , !g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.RunSim(lbm.g);
}
}
