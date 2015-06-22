#include <iostream>
#include <stdexcept>  // runtime_error
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"
#include "UnitTest++.h"

SUITE(TestException)
{
TEST(ConstructorException)
{
//  CHECK_THROW(LatticeBoltzmann lbm, std::runtime_error);
}
}

SUITE(TestFunctionality)
{
  static const double zero_tol = 1e-20;
  static const double loose_tol = 1e-5;
  static const std::size_t g_ny = 6;
  static const std::size_t g_nx = 7;
  static const double g_dx = 0.0316;
  static const double g_dt = 0.001;
  static const double g_t_total = 1.0;
  static const double g_d_coeff = 0.2;
  static const double g_k_visco = 0.2;
  static const double g_rho0_f = 1.1;
  static const double g_rho0_g = 1.2;
  static const std::vector<double> g_u0{1.3, 1.4};
  static const std::vector<std::vector<std::size_t>> g_src_pos_f = {{1, 1}};
  static const std::vector<std::vector<double>> g_src_str_f = {{1.5, 1.6}};
  static const std::vector<std::vector<std::size_t>> g_src_pos_g = {{2, 2}};
  static const std::vector<double> g_src_str_g = {1.7};
  static const std::vector<std::vector<std::size_t>> g_obs_pos;
  static const bool g_is_ns = true;
  static const bool g_is_cd = true;
  static const bool g_is_instant = true;
  static const bool g_no_obstacles = false;

TEST(Constructor)
{
  // formatting function calling like this because input parameters may still
  // change and it's easier to edit this way
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f
    , g_u0);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g
    , g_u0);
  LatticeBoltzmann lbm(g_t_total
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
//  cd.Collide();
  cd.ApplyForce();
  lbm.Print(cd.source_);
  lbm.Print(1, ns.source_);
}

TEST(InitObstacles)
{
  std::vector<bool> obstacle(g_nx * g_ny, false);
  std::vector<std::vector<std::size_t>> obs_pos = {{1, 2}, {3, 4}};
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f
    , g_u0);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g
    , g_u0);
  LatticeBoltzmann lbm(g_t_total
    , obs_pos
    , g_is_ns
    , g_is_cd
    , g_is_instant
    , !g_no_obstacles
    , lm
    , ns
    , cd);
  std::size_t index = 0;
  std::size_t counter = 0;
  for (auto obs : lbm.obstacles) {
    if (counter == obs_pos[index][1] * g_nx + obs_pos[index][0]) {
      CHECK_EQUAL(true, obs);
      if (index < obs_pos.size() - 1) ++index;
    }
    else {
      CHECK_EQUAL(false, obs);
    }
    ++counter;
  }  // obs
}

TEST(InitDensity)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f
    , g_u0);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g
    , g_u0);
  LatticeBoltzmann lbm(g_t_total
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  auto rho_f = lbm.GetRhoF();
  auto rho_g = lbm.GetRhoG();
  for (auto n = 0u; n < rho_f.size(); ++n) {
    CHECK_CLOSE(g_rho0_f, rho_f[n], zero_tol);
    CHECK_CLOSE(g_rho0_g, rho_g[n], zero_tol);
  }  // n
}

TEST(InitVelocity)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f
    , g_u0);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g
    , g_u0);
  LatticeBoltzmann lbm(g_t_total
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  LatticeBoltzmann lbm2(g_t_total
    , g_obs_pos
    , !g_is_ns
    , g_is_cd
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  auto u = lbm.GetVelocity();
  auto u2 = lbm2.GetVelocity();
  for (auto n = 0u; n < u.size(); ++n) {
    CHECK_CLOSE(g_u0[0], u[n][0], zero_tol);
    CHECK_CLOSE(g_u0[1], u[n][1], zero_tol);
    CHECK_CLOSE(g_u0[0], u2[n][0], zero_tol);
    CHECK_CLOSE(g_u0[1], u2[n][1], zero_tol);
  }  // n
}

TEST(ComputeEq)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f
    , g_u0);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g
    , g_u0);
  LatticeBoltzmann lbm(g_t_total
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  std::vector<double> expected_cd = {0.53041,
                                     0.15007, 0.15150, 0.11716, 0.11606,
                                     0.04279, 0.03347, 0.02570, 0.03284};
  std::vector<double> expected_ns = {0.48621,
                                     0.13757, 0.13888, 0.10740, 0.10639,
                                     0.03922, 0.03068, 0.02356, 0.03010};
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    for (auto i = 0; i < 9; ++i) {
      CHECK_CLOSE(expected_ns[i], ns.lattice_eq[n][i], loose_tol);
      CHECK_CLOSE(expected_cd[i], cd.lattice_eq[n][i], loose_tol);
    }  // i
  }  // n
}
}
