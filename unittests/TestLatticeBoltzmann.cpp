#include <iostream>
#include <stdexcept>  // runtime_error
#include <vector>
#include "Algorithm.hpp"
#include "BounceBackNodes.hpp"
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "CollisionNSF.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "Printing.hpp"
#include "StreamPeriodic.hpp"
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
enum Directions {
  E = 1,
  N,
  W,
  S,
  NE,
  NW,
  SW,
  SE
};
static const double zero_tol = 1e-20;
static const double loose_tol = 1e-5;
static const std::size_t g_ny = 6;
static const std::size_t g_nx = 8;
static const double g_dx = 0.0316;
static const double g_dt = 0.001;
static const double g_t_total = 1.0;
static const double g_d_coeff = 0.2;
static const double g_k_visco = 0.2;
static const double g_rho0_f = 1.1;
static const double g_rho0_g = 1.2;
static const std::vector<double> g_u0 = {1.3, 1.4};
static const std::vector<std::vector<std::size_t>> g_src_pos_f = {{1, 1},
    {2, 3}};
static const std::vector<std::vector<double>> g_src_str_f = {{1.5, 1.6},
    {1.7, 1.8}};
static const std::vector<std::vector<std::size_t>> g_src_pos_g = {{2, 2},
    {3, 4}};
static const std::vector<double> g_src_str_g = {1.9, 2.0};
static const std::vector<std::vector<std::size_t>> g_obs_pos;
static const double g_u_lid = 0.01;
static const bool g_is_prestream = true;
static const bool g_is_ns = true;
static const bool g_is_cd = true;
static const bool g_is_taylor = true;
static const bool g_is_lid = true;
static const bool g_is_instant = true;
static const bool g_no_obstacles = false;
/*
TEST(InitObstacles)
{
  std::vector<bool> obstacle(g_nx * g_ny, false);
  std::vector<std::vector<std::size_t>> obs_pos = {{1, 2}, {3, 4}};
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
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
*/
TEST(InitDensity)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    CHECK_CLOSE(g_rho0_f, ns.rho[n], zero_tol);
    CHECK_CLOSE(g_rho0_f, nsf.rho[n], zero_tol);
    CHECK_CLOSE(g_rho0_g, cd.rho[n], zero_tol);
  }  // n
}

TEST(InitVelocity)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  for (auto n = 0u; n < lm.u.size(); ++n) {
    CHECK_CLOSE(g_u0[0], lm.u[n][0], zero_tol);
    CHECK_CLOSE(g_u0[1], lm.u[n][1], zero_tol);
  }  // n
}

TEST(ComputeEq)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  std::vector<double> expected_cd = {0.53041,
                                     0.15007, 0.15150, 0.11716, 0.11606,
                                     0.04279, 0.03347, 0.02570, 0.03284};
  std::vector<double> expected_ns = {0.48621,
                                     0.13757, 0.13888, 0.10740, 0.10639,
                                     0.03922, 0.03068, 0.02356, 0.03010};
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    for (auto i = 0; i < 9; ++i) {
      // can just check the member in cd/ns/nsf since it's using reference
      // instead of creating a new copy
      CHECK_CLOSE(expected_ns[i], ns.edf[n][i], loose_tol);
      CHECK_CLOSE(expected_ns[i], nsf.edf[n][i], loose_tol);
      CHECK_CLOSE(expected_cd[i], cd.edf[n][i], loose_tol);
    }  // i
  }  // n
}

TEST(InitDistributionFunctionLattice)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  LatticeBoltzmann ff(g_obs_pos
    , g_no_obstacles
    , lm
    , nsf);
  LatticeBoltzmann g(g_obs_pos
    , g_no_obstacles
    , lm
    , cd);
  std::vector<double> expected_cd = {0.53041,
                                     0.15007, 0.15150, 0.11716, 0.11606,
                                     0.04279, 0.03347, 0.02570, 0.03284};
  std::vector<double> expected_ns = {0.48621,
                                     0.13757, 0.13888, 0.10740, 0.10639,
                                     0.03922, 0.03068, 0.02356, 0.03010};
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    for (auto i = 0; i < 9; ++i) {
      CHECK_CLOSE(expected_ns[i], f.df[n][i], loose_tol);
      CHECK_CLOSE(expected_ns[i], ff.df[n][i], loose_tol);
      CHECK_CLOSE(expected_cd[i], g.df[n][i], loose_tol);
    }  // i
  }  // n
}

TEST(InitSourceMultiplePosition)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  std::size_t i_ns = 0;
  std::size_t i_cd = 0;
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    if (n == g_src_pos_f[i_ns][1] * g_nx + g_src_pos_f[i_ns][0]) {
      CHECK_CLOSE(g_src_str_f[i_ns][0], nsf.source[n][0], zero_tol);
      CHECK_CLOSE(g_src_str_f[i_ns][1], nsf.source[n][1], zero_tol);
      if (i_ns < g_src_pos_f.size() - 1) ++i_ns;
    }
    else {
      CHECK_CLOSE(0.0, nsf.source[n][0], zero_tol);
      CHECK_CLOSE(0.0, nsf.source[n][1], zero_tol);
    }
    if (n == g_src_pos_g[i_cd][1] * g_nx + g_src_pos_g[i_cd][0]) {
      CHECK_CLOSE(g_src_str_g[i_cd], cd.source[n], zero_tol);
      if (i_cd < g_src_pos_g.size() - 1) ++i_cd;
    }
    else {
      CHECK_CLOSE(0.0, cd.source[n], zero_tol);
    }
  }  // n
}

TEST(CollideWithMixedSource)
{
  // uses global source lattice, get different values at different nodes due to
  // presence of source at that node, nodes without source covers case for
  // CollideNoSource as well
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  LatticeBoltzmann ff(g_obs_pos
    , g_no_obstacles
    , lm
    , nsf);
  LatticeBoltzmann g(g_obs_pos
    , g_no_obstacles
    , lm
    , cd);
  // checks collision result for both cd and ns
  // in expected, index 0 for first source strength, 1 for second source
  // strength (since global pos has 2 positions), 2 for no source (so that it
  // matches source index)
  std::vector<std::vector<double>> expected_cd1 =
      {{0.57428,
        0.22817, 0.22947, 0.19825, 0.19724,
        0.13055, 0.12208, 0.11502, 0.12150},
       {0.57432,
        0.22818, 0.22948, 0.19826, 0.19725,
        0.13056, 0.12208, 0.11502, 0.12151},
       {0.57343,
        0.22795, 0.22924, 0.19805, 0.19705,
        0.13049, 0.12203, 0.11497, 0.12145}};
  std::vector<std::vector<double>> expected_ns1 =
      {{0.53328,
        0.21660, 0.21779, 0.18917, 0.18825,
        0.12726, 0.11949, 0.11302, 0.11896},
       {0.53328,
        0.21660, 0.21779, 0.18917, 0.18825,
        0.12726, 0.11949, 0.11302, 0.11896},
       {0.53328,
        0.21659, 0.21778, 0.18918, 0.18826,
        0.12725, 0.11949, 0.11302, 0.11897}};
  std::vector<std::vector<double>> expected_cd2 =
      {{0.53527,
        0.15745, 0.15887, 0.12479, 0.12369,
        0.05089, 0.04164, 0.03393, 0.04101},
       {0.53532,
        0.15747, 0.15888, 0.12480, 0.12370,
        0.05089, 0.04164, 0.03393, 0.04102},
       {0.53435,
        0.15721, 0.15862, 0.12457, 0.12348,
        0.05083, 0.04158, 0.03388, 0.04095}};
  std::vector<std::vector<double>> expected_ns2 =
      {{0.49052,
        0.14482, 0.14612, 0.11488, 0.11388,
        0.04730, 0.03882, 0.03175, 0.03824},
       {0.49052,
        0.14482, 0.14612, 0.11488, 0.11388,
        0.04730, 0.03882, 0.03175, 0.03824},
       {0.49052,
        0.14481, 0.14611, 0.11489, 0.11389,
        0.04729, 0.03882, 0.03176, 0.03824}};
  // Set distribution function to have different value from equilibrium
  // distribution function so Collide can produce changed result
  f.df.assign(g_nx * g_ny, std::vector<double>(9, 1.0));
  ff.df.assign(g_nx * g_ny, std::vector<double>(9, 1.0));
  g.df.assign(g_nx * g_ny, std::vector<double>(9, 1.0));
  // First collision
  ns.Collide(f.df);
  nsf.Collide(ff.df);
  cd.Collide(g.df);
  for (auto i = 0; i < 9; ++i) {
    std::size_t i_ns = 0;
    std::size_t i_cd = 0;
    for (auto n = 0u; n < g_nx * g_ny; ++n) {
      CHECK_CLOSE(expected_ns1[2][i], f.df[n][i], loose_tol);
      if (n == g_src_pos_f[i_ns][1] * g_nx + g_src_pos_f[i_ns][0]) {
        CHECK_CLOSE(expected_ns1[i_ns][i], ff.df[n][i], loose_tol);
        if (i_ns < g_src_pos_f.size() - 1) ++i_ns;
      }
      else {
        CHECK_CLOSE(expected_ns1[2][i], ff.df[n][i], loose_tol);
      }
      if (n == g_src_pos_g[i_cd][1] * g_nx + g_src_pos_g[i_cd][0]) {
        CHECK_CLOSE(expected_cd1[i_cd][i], g.df[n][i], loose_tol);
        if (i_cd < g_src_pos_g.size() - 1) ++i_cd;
      }
      else {
        CHECK_CLOSE(expected_cd1[2][i], g.df[n][i], loose_tol);
      }
    }  // n
  }  // i
  // Second collision
  ns.Collide(f.df);
  nsf.Collide(ff.df);
  cd.Collide(g.df);
  for (auto i = 0; i < 9; ++i) {
    std::size_t i_ns = 0;
    std::size_t i_cd = 0;
    for (auto n = 0u; n < g_nx * g_ny; ++n) {
      CHECK_CLOSE(expected_ns2[2][i], f.df[n][i], loose_tol);
      if (n == g_src_pos_f[i_ns][1] * g_nx + g_src_pos_f[i_ns][0]) {
        CHECK_CLOSE(expected_ns2[i_ns][i], ff.df[n][i], loose_tol);
        if (i_ns < g_src_pos_f.size() - 1) ++i_ns;
      }
      else {
        CHECK_CLOSE(expected_ns2[2][i], ff.df[n][i], loose_tol);
      }
      if (n == g_src_pos_g[i_cd][1] * g_nx + g_src_pos_g[i_cd][0]) {
        CHECK_CLOSE(expected_cd2[i_cd][i], g.df[n][i], loose_tol);
        if (i_cd < g_src_pos_g.size() - 1) ++i_cd;
      }
      else {
        CHECK_CLOSE(expected_cd2[2][i], g.df[n][i], loose_tol);
      }
    }  // n
  }  // i
}

TEST(ComputeU)
{
  // uses global source lattice, get different values at different nodes due to
  // presence of source at that node, nodes without source covers case for
  // ComputeUNoSource as well
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  LatticeD2Q9 lmf(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  CollisionNSF nsf(lmf
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  LatticeBoltzmann ff(g_obs_pos
    , g_no_obstacles
    , lmf
    , nsf);
  f.df.assign(g_nx * g_ny, {0, 1, 2, 3, 4, 5, 6, 7, 8});
  ff.df.assign(g_nx * g_ny, {0, 1, 2, 3, 4, 5, 6, 7, 8});
  lm.u = ns.ComputeU(f.df);
  lmf.u = nsf.ComputeU(ff.df);
  std::size_t i_ns = 0;
  std::vector<std::vector<double>> expected{{-57.45380, -172.36284},
      {-57.45370, -172.36274}, {-57.45455, -172.36364}};
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    CHECK_CLOSE(expected[2][0], lm.u[n][0], loose_tol);
    CHECK_CLOSE(expected[2][1], lm.u[n][1], loose_tol);
    if (n == g_src_pos_f[i_ns][1] * g_nx + g_src_pos_f[i_ns][0]) {
      CHECK_CLOSE(expected[i_ns][0], lmf.u[n][0], loose_tol);
      CHECK_CLOSE(expected[i_ns][1], lmf.u[n][1], loose_tol);
      if (i_ns < g_src_pos_f.size() - 1) ++i_ns;
    }
    else {
      CHECK_CLOSE(expected[2][0], lmf.u[n][0], loose_tol);
      CHECK_CLOSE(expected[2][1], lmf.u[n][1], loose_tol);
    }
  }  // n
}

TEST(ComputeRho)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  LatticeBoltzmann ff(g_obs_pos
    , g_no_obstacles
    , lm
    , nsf);
  LatticeBoltzmann g(g_obs_pos
    , g_no_obstacles
    , lm
    , cd);
  f.df.assign(g_nx * g_ny, {0, 1, 2, 3, 4, 5, 6, 7, 8});
  ff.df.assign(g_nx * g_ny, {0, 1, 2, 3, 4, 5, 6, 7, 8});
  g.df.assign(g_nx * g_ny, {0, 1, 2, 3, 4, 5, 6, 7, 8});
  ns.rho = ns.ComputeRho(f.df);
  nsf.rho = nsf.ComputeRho(ff.df);
  cd.rho = cd.ComputeRho(g.df);
  auto expected = 36.0;
  for (auto rho : ns.rho) CHECK_CLOSE(expected, rho, loose_tol);
  for (auto rho : nsf.rho) CHECK_CLOSE(expected, rho, loose_tol);
  for (auto rho : cd.rho) CHECK_CLOSE(expected, rho, loose_tol);
}
/*
TEST(BoundaryPeriodic)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.f.assign(g_nx * g_ny, {0, 1, 2, 3, 4, 5, 6, 7, 8});
  lbm.boundary_f = lbm.BoundaryCondition(lbm.f);
  for (auto y = 0u; y < g_ny; ++y) {
    auto n = y * g_nx;
    CHECK_CLOSE(lbm.f[n + g_nx - 1][E], lbm.boundary_f[y][E], zero_tol);
    CHECK_CLOSE(lbm.f[n + g_nx - 1][NE], lbm.boundary_f[y][NE],
        zero_tol);
    CHECK_CLOSE(lbm.f[n + g_nx - 1][SE], lbm.boundary_f[y][SE],
        zero_tol);
    CHECK_CLOSE(lbm.f[n][W], lbm.boundary_f[y + g_ny][W], zero_tol);
    CHECK_CLOSE(lbm.f[n][NW], lbm.boundary_f[y + g_ny][NW], zero_tol);
    CHECK_CLOSE(lbm.f[n][SW], lbm.boundary_f[y + g_ny][SW], zero_tol);
  }  // y
}

TEST(BoundaryPeriodicTaylor)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , g_is_taylor
    , !g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.f.assign(g_nx * g_ny, {0, 1, 2, 3, 4, 5, 6, 7, 8});
  lbm.boundary_f = lbm.BoundaryCondition(lbm.f);
  auto top = 2 * g_ny;
  auto bottom = 2 * g_ny + g_nx;
  for (auto x = 0u; x < g_nx; ++x) {
    auto n = (g_ny - 1) * g_nx;
    CHECK_CLOSE(lbm.f[x][S], lbm.boundary_f[top + x][S], zero_tol);
    CHECK_CLOSE(lbm.f[x][SW], lbm.boundary_f[top + x][SW], zero_tol);
    CHECK_CLOSE(lbm.f[x][SE], lbm.boundary_f[top + x][SE], zero_tol);
    CHECK_CLOSE(lbm.f[n + x][N], lbm.boundary_f[bottom + x][N], zero_tol);
    CHECK_CLOSE(lbm.f[n + x][NW], lbm.boundary_f[bottom + x][NW], zero_tol);
    CHECK_CLOSE(lbm.f[n + x][NE], lbm.boundary_f[bottom + x][NE], zero_tol);
  }  // y
}
*/
TEST(BoundaryBounceback)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  BounceBackNodes bbns(g_is_prestream
    , ns
    , lm);
  BounceBackNodes bbnsf(g_is_prestream
    , nsf
    , lm);
  BounceBackNodes bbcd(g_is_prestream
    , cd
    , lm);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  LatticeBoltzmann ff(g_obs_pos
    , g_no_obstacles
    , lm
    , nsf);
  LatticeBoltzmann g(g_obs_pos
    , g_no_obstacles
    , lm
    , cd);
  std::vector<double> nums = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> bb_nums = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  for (auto &node : f.df) node = nums;
  for (auto &node : ff.df) node = nums;
  for (auto &node : g.df) node = nums;
  for (auto x = 0u; x < g_nx; ++x) {
    bbns.AddNode(x, 0);
    bbnsf.AddNode(x, 1);
    bbcd.AddNode(x, g_ny - 1);
  }
  bbns.UpdateNodes(f.df);
  bbnsf.UpdateNodes(ff.df);
  bbcd.UpdateNodes(g.df);

  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    CHECK_CLOSE(ns.skip[n] ? bb_nums[0] : nums[0], f.df[n][0], zero_tol);
    CHECK_CLOSE(nsf.skip[n] ? bb_nums[0] : nums[0], ff.df[n][0], zero_tol);
    CHECK_CLOSE(cd.skip[n] ? bb_nums[0] : nums[0], g.df[n][0], zero_tol);
  }
}
/*
TEST(BoundaryCorner)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.f.assign(g_nx * g_ny, {0, 1, 2, 3, 4, 5, 6, 7, 8});
  lbm.boundary_f = lbm.BoundaryCondition(lbm.f);
  auto corner = 2 * g_nx + 2 * g_ny;
  CHECK_CLOSE(lbm.f[0][SW], lbm.boundary_f[corner][NE], zero_tol);
  CHECK_CLOSE(lbm.f[g_nx - 1][SE], lbm.boundary_f[corner + 1][NW], zero_tol);
  CHECK_CLOSE(lbm.f[(g_ny - 1) * g_nx][NW], lbm.boundary_f[corner + 2][SE],
      zero_tol);
  CHECK_CLOSE(lbm.f[g_ny * g_nx - 1][NE], lbm.boundary_f[corner + 3][SW],
      zero_tol);
}
*/
TEST(StreamPeriodicHorizontal)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamPeriodic sp(lm);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  std::vector<double> first_three = {1, 2, 1, 0, 1, 2, 0, 0, 2};
  std::vector<double> second_three = {4, 5, 4, 3, 4, 5, 3, 3, 5};
  std::vector<std::vector<double>> nodes = {first_three, second_three};
  std::vector<double> first_result = {1, 5, 1, 3, 1, 5, 3, 3, 5};
  std::vector<double> second_result = {4, 2, 4, 0, 4, 2, 0, 0, 2};
  int counter = 0;
  for (auto &node : f.df) {
    node = nodes[counter++ % 2];
  }  // node
  f.df = sp.Stream(f.df);
  for (auto node : f.df) {
    for (auto i = 0u; i < 9; ++i) {
      if ((node[0] - 1) < zero_tol) {
        CHECK_CLOSE(first_result[i], node[i], zero_tol);
      }
      else {
        CHECK_CLOSE(second_result[i], node[i], zero_tol);
      }
    }  // i
  }  // lat
}

TEST(StreamPeriodicVertical)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamPeriodic sp(lm);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  std::vector<double> first_three = {1, 1, 0, 1, 2, 0, 0, 2, 2};
  std::vector<double> second_three = {4, 4, 3, 4, 5, 3, 3, 5, 5};
  std::vector<std::vector<double>> nodes = {first_three, second_three};
  std::vector<double> first_result = {1, 1, 3, 1, 5, 3, 3, 5, 5};
  std::vector<double> second_result = {4, 4, 0, 4, 2, 0, 0, 2, 2};
  int counter = 0;
  for (auto &node : f.df) node = nodes[(counter++ / g_nx) % 2];
  f.df = sp.Stream(f.df);
  for (auto node : f.df) {
    for (auto i = 0u; i < 9; ++i) {
      if ((node[0] - 1) < zero_tol) {
        CHECK_CLOSE(first_result[i], node[i], zero_tol);
      }
      else {
        CHECK_CLOSE(second_result[i], node[i], zero_tol);
      }
    }  // i
  }  // lat
}

TEST(StreamPeriodicDiagonalNESW)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamPeriodic sp(lm);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  std::vector<double> zeroes(9, 0);
  std::vector<double> ones(9, 1);
  std::vector<double> twos(9, 2);
  std::vector<std::vector<double>> nodes = {zeroes, ones, twos};
  std::vector<double> result_ne = {1, 2, 0};
  std::vector<double> result_sw = {2, 0, 1};
  for (auto y = 0u; y < g_ny; ++y) {
    for (auto x = 0u; x < g_nx; ++x) {
      auto n = y * g_nx + x;
      f.df[n] = nodes[(x % 3 + (y + 2) % 3) % 3];
    }  // x
  }  // y
  f.df = sp.Stream(f.df);
  auto n = 0;
  for (auto node : f.df) {
    CHECK_CLOSE((n % g_nx == 0) ? node[0] : result_ne[node[0]], node[NE],
        zero_tol);
    CHECK_CLOSE((n % g_nx == g_nx - 1) ? node[0] : result_sw[node[0]], node[SW],
        zero_tol);
    ++n;
  }  // n
}

TEST(StreamPeriodicDiagonalNWSE)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamPeriodic sp(lm);
  LatticeBoltzmann f(g_obs_pos
    , g_no_obstacles
    , lm
    , ns);
  std::vector<double> zeroes(9, 0);
  std::vector<double> ones(9, 1);
  std::vector<double> twos(9, 2);
  std::vector<std::vector<double>> nodes = {zeroes, ones, twos};
  std::vector<double> result_nw = {1, 2, 0};
  std::vector<double> result_se = {2, 0, 1};
  for (auto y = 0u; y < g_ny; ++y) {
    for (auto x = 0u; x < g_nx; ++x) {
      auto n = y * g_nx + x;
      f.df[n] = nodes[(2 - x % 3 + (y + 2) % 3) % 3];
    }  // x
  }  // y
  f.df = sp.Stream(f.df);
  auto n = 0;
  for (auto node : f.df) {
    CHECK_CLOSE((n % g_nx == g_nx - 1) ? node[0] : result_nw[node[0]], node[NW],
        zero_tol);
    CHECK_CLOSE((n % g_nx == 0) ? node[0] : result_se[node[0]], node[SE],
        zero_tol);
    ++n;
  }  // n
}
/*
TEST(BoundaryLidLid)
{
  auto u_lid = 0.01;
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  auto value = 0.0;
  for (auto &node : lbm.f) {
    for (auto &i : node) {
      i = value;
      value += 1.0;
    }
  }
  auto f0 = lbm.f;
  lbm.BoundaryLid(lbm.f);
  // check left and right boundary
  for (auto y = 1u; y < g_ny - 1; ++y) {
    auto left = y * g_nx;
    auto right = left + g_nx - 1;
    CHECK_CLOSE(lbm.f[left][W], lbm.f[left][E], zero_tol);
    CHECK_CLOSE(lbm.f[right][E], lbm.f[right][W], zero_tol);
    CHECK_CLOSE(lbm.f[left][SW], lbm.f[left][NE], zero_tol);
    CHECK_CLOSE(lbm.f[left][NW], lbm.f[left][SE], zero_tol);
    CHECK_CLOSE(lbm.f[right][SE], lbm.f[right][NW], zero_tol);
    CHECK_CLOSE(lbm.f[right][NE], lbm.f[right][SW], zero_tol);
  }  // y
  // check top and bottom boundary
  for (auto x = 1u; x < g_nx - 1; ++x) {
    auto top = (g_ny - 1) * g_nx + x;
    auto eq_diff_bottom = 0.5 * (lbm.f[x][E] - lbm.f[x][W]);
    auto eq_diff_top = 0.5 * (lbm.f[top][E] - lbm.f[top][W]);
    auto rho_top = (f0[top][0] + f0[top][E] + f0[top][W] + 2.0 *
        (f0[top][S] + f0[top][SW] + f0[top][SE]));
    CHECK_CLOSE(lbm.f[x][S], lbm.f[x][N], zero_tol);
    CHECK_CLOSE(lbm.f[x][SW], lbm.f[x][NE], zero_tol);
    CHECK_CLOSE(lbm.f[x][SE], lbm.f[x][NW], zero_tol);
    CHECK_CLOSE(lbm.f[top][N], lbm.f[top][S], zero_tol);
    CHECK_CLOSE(lbm.f[top][NW] - eq_diff_top + 0.5 * rho_top * u_lid,
        lbm.f[top][SE], loose_tol);
    CHECK_CLOSE(lbm.f[top][NE] + eq_diff_top - 0.5 * rho_top * u_lid,
        lbm.f[top][SW], loose_tol);
  }  // x
}

TEST(BoundaryCornerZouHeLid)
{
  auto u_lid = 0.01;
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  auto value = 0.0;
  for (auto &node : lbm.f) {
    for (auto &i : node) {
      i = value;
      value += 1.0;
    }  // i
  }  // node
  auto top_left = (g_ny - 1) * g_nx;
  auto top_right = g_ny * g_nx - 1;
  lbm.BoundaryLid(lbm.f);
  CHECK_CLOSE(lbm.f[0][W], lbm.f[0][E], zero_tol);
  CHECK_CLOSE(lbm.f[0][S], lbm.f[0][N], zero_tol);
  CHECK_CLOSE(lbm.f[0][SW], lbm.f[0][NE], zero_tol);

  CHECK_CLOSE(lbm.f[g_nx - 1][E], lbm.f[g_nx - 1][W], zero_tol);
  CHECK_CLOSE(lbm.f[g_nx - 1][S], lbm.f[g_nx - 1][N], zero_tol);
  CHECK_CLOSE(lbm.f[g_nx - 1][SE], lbm.f[g_nx - 1][NW], zero_tol);

  CHECK_CLOSE(318.00667, lbm.f[top_left][E], loose_tol);
  CHECK_CLOSE(317.0, lbm.f[top_left][S], loose_tol);
  CHECK_CLOSE(321.001667, lbm.f[top_left][SE], loose_tol);
  CHECK_CLOSE(0.00083333, lbm.f[top_left][NE], loose_tol);
  CHECK_CLOSE(-0.00083333, lbm.f[top_left][SW], loose_tol);
  CHECK_CLOSE(958.991667, lbm.f[top_left][0], loose_tol);
  CHECK_CLOSE(369.993333, lbm.f[top_right][W], loose_tol);
  CHECK_CLOSE(371.0, lbm.f[top_right][S], loose_tol);
  CHECK_CLOSE(373.998333, lbm.f[top_right][SW], loose_tol);
  CHECK_CLOSE(-0.00083333, lbm.f[top_right][NW], loose_tol);
  CHECK_CLOSE(0.00083333, lbm.f[top_right][SE], loose_tol);
  CHECK_CLOSE(1127.00833, lbm.f[top_right][0], loose_tol);
}

TEST(StreamHorizontalImmersed)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  std::vector<double> first_three = {1, 2, 1, 0, 1, 2, 0, 0, 2};
  std::vector<double> second_three = {4, 5, 4, 3, 4, 5, 3, 3, 5};
  std::vector<std::vector<double>> nodes = {first_three, second_three};
  std::vector<double> first_result = {1, 5, 1, 3, 1, 5, 3, 3, 5};
  std::vector<double> second_result = {4, 2, 4, 0, 4, 2, 0, 0, 2};
  int counter = 0;
  for (auto &node : lbm.f) {
    node = nodes[(counter / g_nx + counter) % 2];
    ++counter;
  }  // node
  lbm.f = lbm.StreamImmersed(lbm.f);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    if ((lbm.f[n][0] - 1) < zero_tol) {
      CHECK_CLOSE((n % g_nx == 0) ? first_three[E] : first_result[E],
          lbm.f[n][E], zero_tol);
      CHECK_CLOSE((n / g_nx == 0) ? first_three[N] : first_result[N],
          lbm.f[n][N], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1) ? first_three[W] : first_result[W],
          lbm.f[n][W], zero_tol);
      CHECK_CLOSE((n / g_nx == g_ny - 1) ? first_three[S] : first_result[S],
          lbm.f[n][S], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == 0) ? first_three[NE] :
          first_result[NE], lbm.f[n][NE], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == 0) ? first_three[NW] :
          first_result[NW], lbm.f[n][NW], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == g_ny - 1) ?
          first_three[SW] : first_result[SW], lbm.f[n][SW], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == g_ny - 1) ?
          first_three[SE] : first_result[SE], lbm.f[n][SE], zero_tol);
    }
    else {
      CHECK_CLOSE((n % g_nx == 0) ? second_three[E] : second_result[E],
          lbm.f[n][E], zero_tol);
      CHECK_CLOSE((n / g_nx == 0) ? second_three[N] : second_result[N],
          lbm.f[n][N], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1) ? second_three[W] : second_result[W],
          lbm.f[n][W], zero_tol);
      CHECK_CLOSE((n / g_nx == g_ny - 1) ? second_three[S] : second_result[S],
          lbm.f[n][S], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == 0) ? second_three[NE] :
          second_result[NE], lbm.f[n][NE], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == 0) ? second_three[NW] :
          second_result[NW], lbm.f[n][NW], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == g_ny - 1) ?
          second_three[SW] : second_result[SW], lbm.f[n][SW], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == g_ny - 1) ?
          second_three[SE] : second_result[SE], lbm.f[n][SE], zero_tol);
    }
  }  // n
}

TEST(StreamVerticalImmersed)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  std::vector<double> first_three = {1, 1, 0, 1, 2, 0, 0, 2, 2};
  std::vector<double> second_three = {4, 4, 3, 4, 5, 3, 3, 5, 5};
  std::vector<std::vector<double>> nodes = {first_three, second_three};
  std::vector<double> first_result = {1, 1, 3, 1, 5, 3, 3, 5, 5};
  std::vector<double> second_result = {4, 4, 0, 4, 2, 0, 0, 2, 2};
  int counter = 0;
  for (auto &node : lbm.f) node = nodes[(counter++ / g_nx) % 2];
  lbm.f = lbm.StreamImmersed(lbm.f);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    if ((lbm.f[n][0] - 1) < zero_tol) {
      CHECK_CLOSE((n % g_nx == 0) ? first_three[E] : first_result[E],
          lbm.f[n][E], zero_tol);
      CHECK_CLOSE((n / g_nx == 0) ? first_three[N] : first_result[N],
          lbm.f[n][N], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1) ? first_three[W] : first_result[W],
          lbm.f[n][W], zero_tol);
      CHECK_CLOSE((n / g_nx == g_ny - 1) ? first_three[S] : first_result[S],
          lbm.f[n][S], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == 0) ? first_three[NE] :
          first_result[NE], lbm.f[n][NE], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == 0) ? first_three[NW] :
          first_result[NW], lbm.f[n][NW], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == g_ny - 1) ?
          first_three[SW] : first_result[SW], lbm.f[n][SW], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == g_ny - 1) ?
          first_three[SE] : first_result[SE], lbm.f[n][SE], zero_tol);
    }
    else {
      CHECK_CLOSE((n % g_nx == 0) ? second_three[E] : second_result[E],
          lbm.f[n][E], zero_tol);
      CHECK_CLOSE((n / g_nx == 0) ? second_three[N] : second_result[N],
          lbm.f[n][N], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1) ? second_three[W] : second_result[W],
          lbm.f[n][W], zero_tol);
      CHECK_CLOSE((n / g_nx == g_ny - 1) ? second_three[S] : second_result[S],
          lbm.f[n][S], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == 0) ? second_three[NE] :
          second_result[NE], lbm.f[n][NE], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == 0) ? second_three[NW] :
          second_result[NW], lbm.f[n][NW], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == g_ny - 1) ?
          second_three[SW] : second_result[SW], lbm.f[n][SW], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == g_ny - 1) ?
          second_three[SE] : second_result[SE], lbm.f[n][SE], zero_tol);
    }
  }  // n
}

TEST(StreamDiagonalNESWImmersed)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  std::vector<double> ones(9, 1);
  std::vector<double> twos(9, 2);
  std::vector<double> threes(9, 3);
  std::vector<std::vector<double>> nodes = {ones, twos, threes};
  std::vector<double> result = {1, 2, 3};
  for (auto y = 0u; y < g_ny; ++y) {
    for (auto x = 0u; x < g_nx; ++x) {
      auto n = y * g_nx + x;
      lbm.f[n] = nodes[(x % 3 + (y + 2) % 3) % 3];
    }  // x
  }  // y
  lbm.f = lbm.StreamImmersed(lbm.f);
  for (auto y = 0u; y < g_ny; ++y) {
    for (auto x = 0u; x < g_nx; ++x) {
      auto n = y * g_nx + x;
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == 0) ? lbm.f[n][0] :
          result[(x % 3 + y % 3) % 3], lbm.f[n][NE], zero_tol);
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == g_ny - 1) ? lbm.f[n][0] :
          result[(x % 3 + y % 3 + 1) % 3], lbm.f[n][SW], zero_tol);
    }  // x
  }  // y
}

TEST(StreamDiagonalNWSEImmerse)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  std::vector<double> ones(9, 1);
  std::vector<double> twos(9, 2);
  std::vector<double> threes(9, 3);
  std::vector<std::vector<double>> nodes = {ones, twos, threes};
  std::vector<double> result = {1, 2, 3};
  for (auto y = 0u; y < g_ny; ++y) {
    for (auto x = 0u; x < g_nx; ++x) {
      auto n = y * g_nx + x;
      lbm.f[n] = nodes[(2 - x % 3 + (y + 2) % 3) % 3];
    }  // x
  }  // y
  lbm.f = lbm.StreamImmersed(lbm.f);
  for (auto y = 0u; y < g_ny; ++y) {
    for (auto x = 0u; x < g_nx; ++x) {
      auto n = y * g_nx + x;
      CHECK_CLOSE((n % g_nx == g_nx - 1 || n / g_nx == 0) ? lbm.f[n][0] :
          result[(2 - x % 3 + y % 3) % 3], lbm.f[n][NW], zero_tol);
      CHECK_CLOSE((n % g_nx == 0 || n / g_nx == g_ny - 1) ? lbm.f[n][0] :
          result[(3 - x % 3 + y % 3) % 3], lbm.f[n][SE], zero_tol);
    }  // x
  }  // y
}

TEST(ObstaclesDoNotReflectIntoObstacles)
{
  std::vector<std::vector<std::size_t>> obs_pos;
  for (auto y = 1u; y < 4u; ++y) {
    obs_pos.push_back({1, y});
    obs_pos.push_back({2, y});
    obs_pos.push_back({3, y});
  }
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , !g_no_obstacles
    , lm
    , ns
    , cd);
  auto obs_center = 2 * g_nx + 2;
  for (auto t = 0; t < 1000; ++t) {
    ns.Collide(lbm.f);
    auto node_f = lbm.f[obs_center];
    lbm.Obstacles(lbm.f);
    // check obstacle center in ns lattice is not changed
    for (auto i = 1u; i < 9u; ++i)
        CHECK_CLOSE(node_f[i], lbm.f[obs_center][i], zero_tol);
    // runs the rest of TakeStep();
    lbm.boundary_f = lbm.BoundaryCondition(lbm.f);
    lbm.f = lbm.Stream(lbm.f, lbm.boundary_f);
    ns.rho = lm.ComputeRho(lbm.f);
    lm.u = lm.ComputeU(lbm.f, ns.rho, ns.source);
    ns.ComputeEq();
    cd.Collide(lbm.g);
    auto node_g = lbm.g[obs_center];
    lbm.Obstacles(lbm.g);
    // check obstacle center in cd lattice is not changed
    for (auto i = 1u; i < 9u; ++i)
        CHECK_CLOSE(node_g[i], lbm.g[obs_center][i], zero_tol);
    lbm.boundary_g = lbm.BoundaryCondition(lbm.g);
    lbm.g = lbm.Stream(lbm.g, lbm.boundary_g);
    cd.rho = lm.ComputeRho(lbm.g);
    cd.ComputeEq();
  }  // t
}

TEST(ObstaclesDoNotReflectLeftRightEdges)
{
  std::vector<std::vector<std::size_t>> obs_pos;
  for (auto y = 0u; y < g_ny; ++y) {
    obs_pos.push_back({0, y});
    obs_pos.push_back({g_nx - 1, y});
  }
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , !g_no_obstacles
    , lm
    , ns
    , cd);
  for (auto t = 0; t < 1000; ++t) {
    ns.Collide(lbm.f);
    std::vector<std::vector<double>> left_edge_unchanged_f;
    std::vector<std::vector<double>> right_edge_unchanged_f;
    for (auto y = 0u; y < g_ny; ++y) {
      auto n = y * g_nx;
      left_edge_unchanged_f.push_back(lbm.f[n]);
      right_edge_unchanged_f.push_back(lbm.f[n + g_nx - 1]);
    }
    lbm.Obstacles(lbm.f);
    // check left right edge obstacle in ns lattice is not changed
    for (auto y = 0u; y < g_ny; ++y) {
      auto n = y * g_nx;
      CHECK_CLOSE(left_edge_unchanged_f[y][W], lbm.f[n][W], zero_tol);
      CHECK_CLOSE(left_edge_unchanged_f[y][NW], lbm.f[n][NW], zero_tol);
      CHECK_CLOSE(left_edge_unchanged_f[y][SW], lbm.f[n][SW], zero_tol);
      CHECK_CLOSE(right_edge_unchanged_f[y][E],
          lbm.f[n + g_nx - 1][E], zero_tol);
      CHECK_CLOSE(right_edge_unchanged_f[y][NE],
          lbm.f[n + g_nx - 1][NE], zero_tol);
      CHECK_CLOSE(right_edge_unchanged_f[y][SE],
          lbm.f[n + g_nx - 1][SE], zero_tol);
    }  // y
    // runs the rest of TakeStep();
    lbm.boundary_f = lbm.BoundaryCondition(lbm.f);
    lbm.f = lbm.Stream(lbm.f, lbm.boundary_f);
    ns.rho = lm.ComputeRho(lbm.f);
    lm.u = lm.ComputeU(lbm.f, ns.rho, ns.source);
    ns.ComputeEq();
    cd.Collide(lbm.g);
    std::vector<std::vector<double>> left_edge_unchanged_g;
    std::vector<std::vector<double>> right_edge_unchanged_g;
    for (auto y = 0u; y < g_ny; ++y) {
      auto n = y * g_nx;
      left_edge_unchanged_g.push_back(lbm.g[n]);
      right_edge_unchanged_g.push_back(lbm.g[n + g_nx - 1]);
    }
    lbm.Obstacles(lbm.g);
    // check left right edge obstacle in cd lattice is not changed
    for (auto y = 0u; y < g_ny; ++y) {
      auto n = y * g_nx;
      CHECK_CLOSE(left_edge_unchanged_g[y][W], lbm.g[n][W], zero_tol);
      CHECK_CLOSE(left_edge_unchanged_g[y][NW], lbm.g[n][NW], zero_tol);
      CHECK_CLOSE(left_edge_unchanged_g[y][SW], lbm.g[n][SW], zero_tol);
      CHECK_CLOSE(right_edge_unchanged_g[y][E],
          lbm.g[n + g_nx - 1][E], zero_tol);
      CHECK_CLOSE(right_edge_unchanged_g[y][NE],
          lbm.g[n + g_nx - 1][NE], zero_tol);
      CHECK_CLOSE(right_edge_unchanged_g[y][SE],
          lbm.g[n + g_nx - 1][SE], zero_tol);
    }
    lbm.boundary_g = lbm.BoundaryCondition(lbm.g);
    lbm.g = lbm.Stream(lbm.g, lbm.boundary_g);
    cd.rho = lm.ComputeRho(lbm.g);
    cd.ComputeEq();
  }  // t
}

TEST(ObstaclesDoNotReflectTopBottomEdges)
{
  std::vector<std::vector<std::size_t>> obs_pos;
  for (auto x = 0u; x < g_nx; ++x) {
    obs_pos.push_back({x, 0});
    obs_pos.push_back({x, g_ny - 1});
  }
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , !g_no_obstacles
    , lm
    , ns
    , cd);
  for (auto t = 0; t < 1000; ++t) {
    ns.Collide(lbm.f);
    std::vector<std::vector<double>> top_edge_unchanged_f;
    std::vector<std::vector<double>> bottom_edge_unchanged_f;
    for (auto x = 0u; x < g_nx; ++x) {
      auto top = (g_ny - 1) * g_nx;
      top_edge_unchanged_f.push_back(lbm.f[top]);
      bottom_edge_unchanged_f.push_back(lbm.f[x]);
    }
    lbm.Obstacles(lbm.f);
    // check top bottom edge obstacle in ns lattice is not changed
    for (auto x = 0u; x < g_nx; ++x) {
      auto top = (g_ny - 1) * g_nx;
      CHECK_CLOSE(top_edge_unchanged_f[x][N], lbm.f[top][N], zero_tol);
      CHECK_CLOSE(top_edge_unchanged_f[x][NW], lbm.f[top][NW], zero_tol);
      CHECK_CLOSE(top_edge_unchanged_f[x][NE], lbm.f[top][NE], zero_tol);
      CHECK_CLOSE(bottom_edge_unchanged_f[x][S], lbm.f[x][S], zero_tol);
      CHECK_CLOSE(bottom_edge_unchanged_f[x][SW], lbm.f[x][SW], zero_tol);
      CHECK_CLOSE(bottom_edge_unchanged_f[x][SE], lbm.f[x][SE], zero_tol);
    }  // y
    // runs the rest of TakeStep();
    lbm.boundary_f = lbm.BoundaryCondition(lbm.f);
    lbm.f = lbm.Stream(lbm.f, lbm.boundary_f);
    ns.rho = lm.ComputeRho(lbm.f);
    lm.u = lm.ComputeU(lbm.f, ns.rho, ns.source);
    ns.ComputeEq();
    cd.Collide(lbm.g);
    std::vector<std::vector<double>> top_edge_unchanged_g;
    std::vector<std::vector<double>> bottom_edge_unchanged_g;
    for (auto x = 0u; x < g_nx; ++x) {
      auto top = (g_ny - 1) * g_nx;
      top_edge_unchanged_g.push_back(lbm.g[top]);
      bottom_edge_unchanged_g.push_back(lbm.g[x]);
    }
    lbm.Obstacles(lbm.g);
    // check top bottom edge obstacle in cd lattice is not changed
    for (auto x = 0u; x < g_nx; ++x) {
      auto top = (g_ny - 1) * g_nx;
      CHECK_CLOSE(top_edge_unchanged_g[x][N], lbm.g[top][N], zero_tol);
      CHECK_CLOSE(top_edge_unchanged_g[x][NW], lbm.g[top][NW], zero_tol);
      CHECK_CLOSE(top_edge_unchanged_g[x][NE], lbm.g[top][NE], zero_tol);
      CHECK_CLOSE(bottom_edge_unchanged_g[x][S], lbm.g[x][S], zero_tol);
      CHECK_CLOSE(bottom_edge_unchanged_g[x][SW], lbm.g[x][SW], zero_tol);
      CHECK_CLOSE(bottom_edge_unchanged_g[x][SE], lbm.g[x][SE], zero_tol);
    }
    lbm.boundary_g = lbm.BoundaryCondition(lbm.g);
    lbm.g = lbm.Stream(lbm.g, lbm.boundary_g);
    cd.rho = lm.ComputeRho(lbm.g);
    cd.ComputeEq();
  }  // t
}

TEST(ObstaclesDoNotReflectCorners)
{
  std::vector<std::vector<std::size_t>> obs_pos = {{0, 0}, {0, g_ny - 1},
      {g_nx - 1, 0}, {g_nx - 1, g_ny - 1}};
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , !g_no_obstacles
    , lm
    , ns
    , cd);
  auto bot_right = g_nx - 1;
  auto top_left = (g_ny - 1) * g_nx;
  auto top_right = (g_ny - 1) * g_nx + g_nx - 1;
  for (auto t = 0; t < 1000; ++t) {
    ns.Collide(lbm.f);
    std::vector<std::vector<double>> corners_unchanged_f = {lbm.f[0],
        lbm.f[bot_right], lbm.f[top_left], lbm.f[top_right]};
    lbm.Obstacles(lbm.f);
    // check corner obstacles in ns lattice is not changed
    CHECK_CLOSE(corners_unchanged_f[0][SW], lbm.f[0][SW], zero_tol);
    CHECK_CLOSE(corners_unchanged_f[1][SE], lbm.f[bot_right][SE], zero_tol);
    CHECK_CLOSE(corners_unchanged_f[2][NW], lbm.f[top_left][NW], zero_tol);
    CHECK_CLOSE(corners_unchanged_f[3][NE], lbm.f[top_right][NE], zero_tol);
    // runs the rest of TakeStep();
    lbm.boundary_f = lbm.BoundaryCondition(lbm.f);
    lbm.f = lbm.Stream(lbm.f, lbm.boundary_f);
    ns.rho = lm.ComputeRho(lbm.f);
    lm.u = lm.ComputeU(lbm.f, ns.rho, ns.source);
    ns.ComputeEq();
    cd.Collide(lbm.g);
    std::vector<std::vector<double>> top_edge_unchanged_g;
    std::vector<std::vector<double>> bottom_edge_unchanged_g;
    std::vector<std::vector<double>> corners_unchanged_g = {lbm.g[0],
        lbm.g[bot_right], lbm.g[top_left], lbm.g[top_right]};
    lbm.Obstacles(lbm.g);
    // check corner obstacles in cd lattice is not changed
    CHECK_CLOSE(corners_unchanged_g[0][SW], lbm.g[0][SW], zero_tol);
    CHECK_CLOSE(corners_unchanged_g[1][SE], lbm.g[bot_right][SE], zero_tol);
    CHECK_CLOSE(corners_unchanged_g[2][NW], lbm.g[top_left][NW], zero_tol);
    CHECK_CLOSE(corners_unchanged_g[3][NE], lbm.g[top_right][NE], zero_tol);
    lbm.boundary_g = lbm.BoundaryCondition(lbm.g);
    lbm.g = lbm.Stream(lbm.g, lbm.boundary_g);
    cd.rho = lm.ComputeRho(lbm.g);
    cd.ComputeEq();
  }  // t
}

TEST(ObstaclesNormalReflect)
{
  std::vector<std::vector<std::size_t>> obs_pos = {{2, 4}};
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , !g_no_obstacles
    , lm
    , ns
    , cd);
  auto n = 4 * g_nx + 2;
  for (auto t = 0; t < 1; ++t) {
    ns.Collide(lbm.f);
    lbm.Obstacles(lbm.f);
    // check directions around the obstacle in ns lattice are reflected
    CHECK_CLOSE(lbm.f[n + g_nx][S], lbm.f[n][N], zero_tol);
    CHECK_CLOSE(lbm.f[n + 1][W], lbm.f[n][E], zero_tol);
    CHECK_CLOSE(lbm.f[n - 1][E], lbm.f[n][W], zero_tol);
    CHECK_CLOSE(lbm.f[n - g_nx][N], lbm.f[n][S], zero_tol);
    CHECK_CLOSE(lbm.f[n + g_nx + 1][SW], lbm.f[n][NE], zero_tol);
    CHECK_CLOSE(lbm.f[n + g_nx - 1][SE], lbm.f[n][NW], zero_tol);
    CHECK_CLOSE(lbm.f[n - g_nx - 1][NE], lbm.f[n][SW], zero_tol);
    CHECK_CLOSE(lbm.f[n - g_nx + 1][NW], lbm.f[n][SE], zero_tol);
    // runs the rest of TakeStep();
    lbm.boundary_f = lbm.BoundaryCondition(lbm.f);
    lbm.f = lbm.Stream(lbm.f, lbm.boundary_f);
    ns.rho = lm.ComputeRho(lbm.f);
    lm.u = lm.ComputeU(lbm.f, ns.rho, ns.source);
    ns.ComputeEq();
    cd.Collide(lbm.g);
    lbm.Obstacles(lbm.g);
    // check directions around the obstacle in cd lattice are reflected
    CHECK_CLOSE(lbm.g[n + g_nx][S], lbm.g[n][N], zero_tol);
    CHECK_CLOSE(lbm.g[n + 1][W], lbm.g[n][E], zero_tol);
    CHECK_CLOSE(lbm.g[n - 1][E], lbm.g[n][W], zero_tol);
    CHECK_CLOSE(lbm.g[n - g_nx][N], lbm.g[n][S], zero_tol);
    CHECK_CLOSE(lbm.g[n + g_nx + 1][SW], lbm.g[n][NE], zero_tol);
    CHECK_CLOSE(lbm.g[n + g_nx - 1][SE], lbm.g[n][NW], zero_tol);
    CHECK_CLOSE(lbm.g[n - g_nx - 1][NE], lbm.g[n][SW], zero_tol);
    CHECK_CLOSE(lbm.g[n - g_nx + 1][NW], lbm.g[n][SE], zero_tol);
    lbm.boundary_g = lbm.BoundaryCondition(lbm.g);
    lbm.g = lbm.Stream(lbm.g, lbm.boundary_g);
    cd.rho = lm.ComputeRho(lbm.g);
    cd.ComputeEq();
  }  // t
}

TEST(InstantSourceToggle)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  LatticeBoltzmann lbm(g_t_total
    , g_u_lid
    , g_obs_pos
    , g_is_ns
    , g_is_cd
    , !g_is_taylor
    , !g_is_lid
    , g_is_instant
    , g_no_obstacles
    , lm
    , ns
    , cd);
  lbm.RunSim();
  // check cd source is zero
  for (auto node : cd.source) CHECK_CLOSE(0.0, node, zero_tol);
  // check ns source is unchanged
  std::size_t i_ns = 0;
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    if (n == g_src_pos_f[i_ns][1] * g_nx + g_src_pos_f[i_ns][0]) {
      CHECK_CLOSE(g_src_str_f[i_ns][0], ns.source[n][0], zero_tol);
      CHECK_CLOSE(g_src_str_f[i_ns][1], ns.source[n][1], zero_tol);
      if (i_ns < g_src_pos_f.size() - 1) ++i_ns;
    }
    else {
      CHECK_CLOSE(0.0, ns.source[n][0], zero_tol);
      CHECK_CLOSE(0.0, ns.source[n][1], zero_tol);
    }
  }  // n
}
*/
}
