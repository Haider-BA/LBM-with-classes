#include <cmath>
#include <iostream>
#include <stdexcept>  // runtime_error
#include <vector>
#include "Algorithm.hpp"
#include "BoundaryNodes.hpp"
#include "BouncebackNodes.hpp"
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "CollisionNSF.hpp"
#include "ImmersedBoundaryMethod.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "Particle.hpp"
#include "ParticleRigid.hpp"
#include "Printing.hpp"
#include "StreamD2Q9.hpp"
#include "StreamPeriodic.hpp"
#include "UnitTest++.h"
#include "ZouHeNodes.hpp"
#include "ZouHePressureNodes.hpp"

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
static const double g_d_coeff = 0.2;
static const double g_k_visco = 0.2;
static const double g_rho0_f = 1.1;
static const double g_rho0_g = 1.2;
static const double g_pi = 3.14159265;
static const std::vector<double> g_u0 = {1.3, 1.4};
static const std::vector<std::vector<std::size_t>> g_src_pos_f = {{1, 1},
    {2, 3}};
static const std::vector<std::vector<double>> g_src_str_f = {{1.5, 1.6},
    {1.7, 1.8}};
static const std::vector<std::vector<std::size_t>> g_src_pos_g = {{2, 2},
    {3, 4}};
static const std::vector<double> g_src_str_g = {1.9, 2.0};
static const bool g_is_prestream = true;
static const bool g_is_modify_stream = true;
static const bool g_is_instant = true;

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
    , g_rho0_g
    , !g_is_instant);
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
    , g_rho0_g
    , !g_is_instant);
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
  StreamPeriodic sp(lm);
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
    , g_rho0_g
    , !g_is_instant);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  LatticeBoltzmann ff(lm
    , nsf
    , sp);
  LatticeBoltzmann g(lm
    , cd
    , sp);
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
    , g_rho0_g
    , !g_is_instant);
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
  StreamPeriodic sp(lm);
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
    , g_rho0_g
    , !g_is_instant);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  LatticeBoltzmann ff(lm
    , nsf
    , sp);
  LatticeBoltzmann g(lm
    , cd
    , sp);
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
  StreamPeriodic sp(lm);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  CollisionNSF nsf(lmf
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  LatticeBoltzmann ff(lm
    , nsf
    , sp);
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
  StreamPeriodic sp(lm);
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
    , g_rho0_g
    , !g_is_instant);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  LatticeBoltzmann ff(lm
    , nsf
    , sp);
  LatticeBoltzmann g(lm
    , cd
    , sp);
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

TEST(BoundaryBounceback)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  StreamPeriodic sp(lm);
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
    , g_rho0_g
    , !g_is_instant);
  BouncebackNodes bbns(lm
    , &ns);
  BouncebackNodes bbnsf(lm
    , &nsf);
  BouncebackNodes bbcd(lm
    , &cd);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  LatticeBoltzmann ff(lm
    , nsf
    , sp);
  LatticeBoltzmann g(lm
    , cd
    , sp);
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
  bbns.UpdateNodes(f.df
    , !g_is_modify_stream);
  bbnsf.UpdateNodes(ff.df
    , !g_is_modify_stream);
  bbcd.UpdateNodes(g.df
    , !g_is_modify_stream);

  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    CHECK_CLOSE(ns.skip[n] ? bb_nums[0] : nums[0], f.df[n][0], zero_tol);
    CHECK_CLOSE(nsf.skip[n] ? bb_nums[0] : nums[0], ff.df[n][0], zero_tol);
    CHECK_CLOSE(cd.skip[n] ? bb_nums[0] : nums[0], g.df[n][0], zero_tol);
  }
}

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
  LatticeBoltzmann f(lm
    , ns
    , sp);
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
  LatticeBoltzmann f(lm
    , ns
    , sp);
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
  LatticeBoltzmann f(lm
    , ns
    , sp);
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
  LatticeBoltzmann f(lm
    , ns
    , sp);
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

TEST(BoundaryZouHeSide)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  ZouHeNodes zhns(lm
    , ns);
  StreamPeriodic sp(lm);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  auto u_lid = 0.1;
  auto v_lid = 0.1;
  std::vector<double> nums = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> right_ans = {1,
                                   2, 3, 1.909569, 5,
                                   6, 10.045216, 4.909569, 9};
  std::vector<double> top_ans = {1,
                                 2, 3, 4, 2.917981,
                                 6, 7, 4.917981, 8.041009};
  std::vector<double> left_ans = {1,
                                  4.099471, 3, 4, 5,
                                  9.099471, 7, 8, 5.950265};
  std::vector<double> bottom_ans = {1,
                                    2, 5.107937, 4, 5,
                                    9.107937, 7.946032, 8, 9};
  for (auto y = 1u; y < g_ny - 1; ++y) {
    auto left = y * g_nx;
    auto right = left + g_nx - 1;
    zhns.AddNode(g_nx - 1, y, u_lid, v_lid);
    zhns.AddNode(0, y, u_lid, v_lid);
    f.df[right] = nums;
    f.df[left] = nums;
  }  // y
  for (auto x = 1u; x < g_nx - 1; ++x) {
    auto top = (g_ny - 1) * g_nx + x;
    zhns.AddNode(x, g_ny - 1, u_lid, v_lid);
    zhns.AddNode(x, 0, u_lid, v_lid);
    f.df[top] = nums;
    f.df[x] = nums;
  }  // x
  zhns.UpdateNodes(f.df, !g_is_modify_stream);
  for (auto i = 0; i < 9; ++i) {
    for (auto y = 1u; y < g_ny - 1; ++y) {
      auto left = y * g_nx;
      auto right = left + g_nx - 1;
      CHECK_CLOSE(right_ans[i], f.df[right][i], loose_tol);
      CHECK_CLOSE(left_ans[i], f.df[left][i], loose_tol);
    }  // y
    for (auto x = 1u; x < g_nx - 1; ++x) {
      auto top = (g_ny - 1) * g_nx + x;
      CHECK_CLOSE(top_ans[i], f.df[top][i], loose_tol);
      CHECK_CLOSE(bottom_ans[i], f.df[x][i], loose_tol);
    }  // x
  }  // i
}

TEST(BoundaryZouHeCorner)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  ZouHeNodes zhns(lm
    , ns);
  StreamPeriodic sp(lm);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  auto u_lid = 0.1;
  auto v_lid = 0.1;
  std::vector<double> nums = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> bottom_left = {10.762658,
                                     4.094937, 5.094937, 4, 5,
                                     8.047468, 0, 8, 0};
  std::vector<double> bottom_right = {13,
                                      2, 5.094937, 1.905063, 5,
                                      0.023734, 9, -0.023734, 9};
  std::vector<double> top_left = {17,
                                  4.094937, 3, 4, 2.905063,
                                  0.023734, 7, -0.023734, 7};
  std::vector<double> top_right = {23.23734,
                                   2, 3, 1.905063, 2.905063,
                                   6, 0, 5.952532, 0};
  for (auto &node : f.df) node = nums;
  ns.rho = ns.ComputeRho(f.df);
  zhns.AddNode(0, 0, u_lid, v_lid);
  zhns.AddNode(g_nx - 1, 0, u_lid, v_lid);
  zhns.AddNode(0, g_ny - 1, u_lid, v_lid);
  zhns.AddNode(g_nx - 1, g_ny - 1, u_lid, v_lid);
  zhns.UpdateNodes(f.df, !g_is_prestream);
  for (auto i = 0; i < 9; ++i) {
    CHECK_CLOSE(bottom_left[i], f.df[0][i], loose_tol);
    CHECK_CLOSE(bottom_right[i], f.df[g_nx - 1][i], loose_tol);
    CHECK_CLOSE(top_left[i], f.df[(g_ny - 1) * g_nx][i], loose_tol);
    CHECK_CLOSE(top_right[i], f.df[g_ny * g_nx - 1][i], loose_tol);
  }
}

TEST(BoundaryZouHePressureSide)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  ZouHePressureNodes zhpns(lm
    , ns);
  StreamPeriodic sp(lm);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  auto rho_node = 55.5;
  std::vector<double> nums = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> right_ans = {1,
                                   2, 3, 10.333333, 5,
                                   6, 12.083333, 7.083333, 9};
  std::vector<double> top_ans = {1,
                                 2, 3, 4, 14,
                                 6, 7, 7.75, 10.75};
  std::vector<double> left_ans = {1,
                                  9.666667, 3, 4, 5,
                                  10.416667, 7, 8, 7.416667};
  std::vector<double> bottom_ans = {1,
                                    2, 8, 4, 5,
                                    9.75, 8.75, 8, 9};
  for (auto y = 1u; y < g_ny - 1; ++y) {
    auto left = y * g_nx;
    auto right = left + g_nx - 1;
    zhpns.AddNode(g_nx - 1, y, rho_node);
    zhpns.AddNode(0, y, rho_node);
    f.df[right] = nums;
    f.df[left] = nums;
  }  // y
  for (auto x = 1u; x < g_nx - 1; ++x) {
    auto top = (g_ny - 1) * g_nx + x;
    zhpns.AddNode(x, g_ny - 1, rho_node);
    zhpns.AddNode(x, 0, rho_node);
    f.df[top] = nums;
    f.df[x] = nums;
  }  // x
  zhpns.UpdateNodes(f.df, !g_is_modify_stream);
  for (auto i = 0; i < 9; ++i) {
    for (auto y = 1u; y < g_ny - 1; ++y) {
      auto left = y * g_nx;
      auto right = left + g_nx - 1;
      CHECK_CLOSE(right_ans[i], f.df[right][i], loose_tol);
      CHECK_CLOSE(left_ans[i], f.df[left][i], loose_tol);
    }  // y
    for (auto x = 1u; x < g_nx - 1; ++x) {
      auto top = (g_ny - 1) * g_nx + x;
      CHECK_CLOSE(top_ans[i], f.df[top][i], loose_tol);
      CHECK_CLOSE(bottom_ans[i], f.df[x][i], loose_tol);
    }  // x
  }  // i
}

TEST(BoundaryZouHePressureCorner)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  ZouHePressureNodes zhpns(lm
    , ns);
  StreamPeriodic sp(lm);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  auto rho_node = 55.5;
  std::vector<double> nums = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> bottom_left = {1,
                                     4, 5, 4, 5,
                                     8, 10.25, 8, 10.25};
  std::vector<double> bottom_right = {1,
                                      2, 5, 2, 5,
                                      11.25, 9, 11.25, 9};
  std::vector<double> top_left = {1,
                                  4, 3, 4, 3,
                                  13.25, 7, 13.25, 7};
  std::vector<double> top_right = {1,
                                   2, 3, 2, 3,
                                   6, 16.25, 6, 16.25};
  for (auto &node : f.df) node = nums;
  ns.rho = ns.ComputeRho(f.df);
  zhpns.AddNode(0, 0, rho_node);
  zhpns.AddNode(g_nx - 1, 0, rho_node);
  zhpns.AddNode(0, g_ny - 1, rho_node);
  zhpns.AddNode(g_nx - 1, g_ny - 1, rho_node);
  zhpns.UpdateNodes(f.df, !g_is_modify_stream);
  for (auto i = 0; i < 9; ++i) {
    CHECK_CLOSE(bottom_left[i], f.df[0][i], loose_tol);
    CHECK_CLOSE(bottom_right[i], f.df[g_nx - 1][i], loose_tol);
    CHECK_CLOSE(top_left[i], f.df[(g_ny - 1) * g_nx][i], loose_tol);
    CHECK_CLOSE(top_right[i], f.df[g_ny * g_nx - 1][i], loose_tol);
  }
}

TEST(HalfwayBouncebackCopyDF)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamD2Q9 sd(lm);
  StreamPeriodic sp(lm);
  BouncebackNodes hwbb(lm
    , &sd);
  BouncebackNodes hwbbsp(lm
    , &sp);
  LatticeBoltzmann f(lm
    , ns
    , sd);
  LatticeBoltzmann ff(lm
    , ns
    , sp);
  std::vector<bool> bounce_back(g_nx * g_ny, false);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    for (auto i = 0u; i < 9; ++i) {
      f.df[n][i] = static_cast<double>(n) + static_cast<double>(i) * 0.1;
      ff.df[n][i] = static_cast<double>(n) + static_cast<double>(i) * 0.1;
    }  // i
  }  // n
  for (auto y = 0u; y < g_ny; ++y) {
    hwbb.AddNode(0, y);
    hwbb.AddNode(g_nx - 1, y);
    hwbbsp.AddNode(0, y);
    hwbbsp.AddNode(g_nx - 1, y);
    bounce_back[y * g_nx] = true;
    bounce_back[y * g_nx + g_nx - 1] = true;
  }
  hwbb.UpdateNodes(f.df
    , !g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , !g_is_modify_stream);
  for (auto node : hwbb.nodes) {
    auto ind = 0;
    for (auto i : node.df_node) CHECK_CLOSE(f.df[node.n][ind++], i, zero_tol);
  }
  for (auto node : hwbbsp.nodes) {
    auto ind = 0;
    for (auto i : node.df_node) CHECK_CLOSE(ff.df[node.n][ind++], i, zero_tol);
  }
}

TEST(StreamHorizontalHalfwayBounceback)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamD2Q9 sd(lm);
  StreamPeriodic sp(lm);
  BouncebackNodes hwbb(lm
    , &sd);
  BouncebackNodes hwbbsp(lm
    , &sp);
  LatticeBoltzmann f(lm
    , ns
    , sd);
  LatticeBoltzmann ff(lm
    , ns
    , sp);
  std::vector<bool> bounce_back(g_nx * g_ny, false);
  for (auto y = 0u; y < g_ny; ++y) {
    hwbb.AddNode(0, y);
    hwbb.AddNode(g_nx - 1, y);
    hwbbsp.AddNode(0, y);
    hwbbsp.AddNode(g_nx - 1, y);
    bounce_back[y * g_nx] = true;
    bounce_back[y * g_nx + g_nx - 1] = true;
  }
  std::vector<double> first_three = {1, 2, 1, 0, 1, 2, 0, 0, 2};
  std::vector<double> second_three = {4, 5, 4, 3, 4, 5, 3, 3, 5};
  std::vector<std::vector<double>> nodes = {first_three, second_three};
  std::vector<double> first_result = {1, 5, 1, 3, 1, 5, 3, 3, 5};
  std::vector<double> second_result = {4, 2, 4, 0, 4, 2, 0, 0, 2};
  std::vector<double> first_top_result = {1, 5, 1, 3, 1, 5, 3, 0, 2};
  std::vector<double> second_top_result = {4, 2, 4, 0, 4, 2, 0, 3, 5};
  std::vector<double> first_bottom_result = {1, 5, 1, 3, 1, 2, 0, 3, 5};
  std::vector<double> second_bottom_result = {4, 2, 4, 0, 4, 5, 3, 0, 2};
  std::vector<double> left_result = {1, 0, 1, 3, 1, 0, 3, 3, 0};
  std::vector<double> right_result = {4, 2, 4, 5, 4, 2, 5, 5, 2};
  std::vector<double> bottom_left_result = {1, 0, 1, 3, 1, 0, 2, 3, 0};
  std::vector<double> bottom_right_result = {4, 2, 4, 5, 4, 3, 5, 5, 2};
  std::vector<double> top_left_result = {1, 0, 1, 3, 1, 0, 3, 2, 0};
  std::vector<double> top_right_result = {4, 2, 4, 5, 4, 2, 5, 5, 3};
  int counter = 0;
  for (auto &node : f.df) node = nodes[counter++ % 2];
  counter = 0;
  for (auto &node : ff.df) node = nodes[counter++ % 2];
  hwbb.UpdateNodes(f.df
    , !g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , !g_is_modify_stream);
  f.df = sd.Stream(f.df);
  ff.df = sp.Stream(ff.df);
  hwbb.UpdateNodes(f.df
    , g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , g_is_modify_stream);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    auto bottom = n / g_nx == 0;
    auto top = n / g_nx == g_ny - 1;
    for (auto i = 0u; i < 9; ++i) {
      if (bounce_back[n]) {
        if (top) {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? top_left_result[i] :
              top_right_result[i], f.df[n][i], zero_tol);
          CHECK_CLOSE(ff.df[n][0] - 1.0 < zero_tol ? top_left_result[i] :
              top_right_result[i], ff.df[n][i], zero_tol);
        }
        else if (bottom) {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? bottom_left_result[i] :
              bottom_right_result[i], f.df[n][i], zero_tol);
          CHECK_CLOSE(ff.df[n][0] - 1.0 < zero_tol ? bottom_left_result[i] :
              bottom_right_result[i], ff.df[n][i], zero_tol);
        }
        else {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? left_result[i] :
              right_result[i], f.df[n][i], zero_tol);
          CHECK_CLOSE(ff.df[n][0] - 1.0 < zero_tol ? left_result[i] :
              right_result[i], ff.df[n][i], zero_tol);
        }
      }
      else {
        if (top) {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? first_top_result[i] :
              second_top_result[i], f.df[n][i], zero_tol);
        }
        else if (bottom) {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? first_bottom_result[i] :
              second_bottom_result[i], f.df[n][i], zero_tol);
        }
        else {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? first_result[i] :
              second_result[i], f.df[n][i], zero_tol);
        }
        CHECK_CLOSE(ff.df[n][0] - 1.0 < zero_tol ? first_result[i] :
            second_result[i], ff.df[n][i], zero_tol);
      }
    }  // i
  }  // n
}

TEST(StreamVerticalHalfwayBounceback)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamD2Q9 sd(lm);
  StreamPeriodic sp(lm);
  BouncebackNodes hwbb(lm
    , &sd);
  BouncebackNodes hwbbsp(lm
    , &sp);
  LatticeBoltzmann f(lm
    , ns
    , sd);
  LatticeBoltzmann ff(lm
    , ns
    , sp);
  std::vector<bool> bounce_back(g_nx * g_ny, false);
  for (auto x = 0u; x < g_nx; ++x) {
    hwbb.AddNode(x, 0);
    hwbb.AddNode(x, g_ny - 1);
    hwbbsp.AddNode(x, 0);
    hwbbsp.AddNode(x, g_ny - 1);
    bounce_back[x] = true;
    bounce_back[(g_ny - 1) * g_nx + x] = true;
  }
  std::vector<double> first_three = {1, 1, 0, 1, 2, 0, 0, 2, 2};
  std::vector<double> second_three = {4, 4, 3, 4, 5, 3, 3, 5, 5};
  std::vector<std::vector<double>> nodes = {first_three, second_three};
  std::vector<double> first_result = {1, 1, 3, 1, 5, 3, 3, 5, 5};
  std::vector<double> second_result = {4, 4, 0, 4, 2, 0, 0, 2, 2};
  std::vector<double> first_left_result = {1, 1, 3, 1, 5, 0, 3, 5, 2};
  std::vector<double> second_left_result = {4, 4, 0, 4, 2, 3, 0, 2, 5};
  std::vector<double> first_right_result = {1, 1, 3, 1, 5, 3, 0, 2, 5};
  std::vector<double> second_right_result = {4, 4, 0, 4, 2, 0, 3, 5, 2};
  std::vector<double> top_result = {4, 4, 0, 4, 3, 0, 0, 3, 3};
  std::vector<double> bottom_result = {1, 1, 2, 1, 5, 2, 2, 5, 5};
  std::vector<double> top_left_result = {4, 4, 0, 4, 3, 5, 0, 3, 3};
  std::vector<double> bottom_left_result = {1, 1, 2, 1, 5, 2, 2, 5, 0};
  std::vector<double> top_right_result = {4, 4, 0, 4, 3, 0, 5, 3, 3};
  std::vector<double> bottom_right_result = {1, 1, 2, 1, 5, 2, 2, 0, 5};
  int counter = 0;
  for (auto &node : f.df) node = nodes[(counter++ / g_nx) % 2];
  counter = 0;
  for (auto &node : ff.df) node = nodes[(counter++ / g_nx) % 2];
  hwbb.UpdateNodes(f.df
    , !g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , !g_is_modify_stream);
  f.df = sd.Stream(f.df);
  ff.df = sp.Stream(ff.df);
  hwbb.UpdateNodes(f.df
    , g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , g_is_modify_stream);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    auto left = n % g_nx == 0;
    auto right = n % g_nx == g_nx - 1;
    for (auto i = 0u; i < 9; ++i) {
      if (bounce_back[n]) {
        if (left) {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? bottom_left_result[i] :
              top_left_result[i], f.df[n][i], zero_tol);
          CHECK_CLOSE(ff.df[n][0] - 1.0 < zero_tol ? bottom_left_result[i] :
              top_left_result[i], ff.df[n][i], zero_tol);
        }
        else if (right) {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? bottom_right_result[i] :
              top_right_result[i], f.df[n][i], zero_tol);
          CHECK_CLOSE(ff.df[n][0] - 1.0 < zero_tol ? bottom_right_result[i] :
              top_right_result[i], ff.df[n][i], zero_tol);
        }
        else {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? bottom_result[i] :
              top_result[i], f.df[n][i], zero_tol);
          CHECK_CLOSE(ff.df[n][0] - 1.0 < zero_tol ? bottom_result[i] :
              top_result[i], ff.df[n][i], zero_tol);
        }
      }
      else {
        if (left) {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? first_left_result[i] :
              second_left_result[i], f.df[n][i], zero_tol);
        }
        else if (right) {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? first_right_result[i] :
              second_right_result[i], f.df[n][i], zero_tol);
        }
        else {
          CHECK_CLOSE(f.df[n][0] - 1.0 < zero_tol ? first_result[i] :
              second_result[i], f.df[n][i], zero_tol);
        }
        CHECK_CLOSE(ff.df[n][0] - 1.0 < zero_tol ? first_result[i] :
              second_result[i], ff.df[n][i], zero_tol);
      }
    }  // i
  }  // n
}

TEST(StreamDiagonalNESWOnGridBounceback)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamD2Q9 sd(lm);
  StreamD2Q9 sd2(lm);
  StreamPeriodic sp(lm);
  StreamPeriodic sp2(lm);
  BouncebackNodes hwbb(lm
    , &sd);
  BouncebackNodes hwbb2(lm
    , &sd2);
  BouncebackNodes hwbbsp(lm
    , &sp);
  BouncebackNodes hwbbsp2(lm
    , &sp2);
  LatticeBoltzmann f(lm
    , ns
    , sd);
  LatticeBoltzmann f2(lm
    , ns
    , sd2);
  LatticeBoltzmann ff(lm
    , ns
    , sp);
  LatticeBoltzmann ff2(lm
    , ns
    , sp2);
  // same results for half-way bounceback with StreamD2Q9 on top & bottom and
  // left & right
  std::vector<bool> bounce_back(g_nx * g_ny, false);
  for (auto y = 0u; y < g_ny; ++y) {
    hwbb.AddNode(0, y);
    hwbb.AddNode(g_nx - 1, y);
    hwbbsp.AddNode(0, y);
    hwbbsp.AddNode(g_nx - 1, y);
    bounce_back[y * g_nx] = true;
    bounce_back[y * g_nx + g_nx - 1] = true;
  }
  std::vector<bool> bounce_back2(g_nx * g_ny, false);
  for (auto x = 0u; x < g_nx; ++x) {
    hwbb2.AddNode(x, 0);
    hwbb2.AddNode(x, g_ny - 1);
    hwbbsp2.AddNode(x, 0);
    hwbbsp2.AddNode(x, g_ny - 1);
    bounce_back2[x] = true;
    bounce_back2[(g_ny - 1) * g_nx + x] = true;
  }
  auto width = g_nx - 1;
  auto height = (g_ny - 1) * g_nx;
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
      f2.df[n] = nodes[(x % 3 + (y + 2) % 3) % 3];
      ff.df[n] = nodes[(x % 3 + (y + 2) % 3) % 3];
      ff2.df[n] = nodes[(x % 3 + (y + 2) % 3) % 3];
    }  // x
  }  // y
  hwbb.UpdateNodes(f.df
    , !g_is_modify_stream);
  hwbb2.UpdateNodes(f2.df
    , !g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , !g_is_modify_stream);
  hwbbsp2.UpdateNodes(ff2.df
    , !g_is_modify_stream);
  f.df = sd.Stream(f.df);
  f2.df = sd2.Stream(f2.df);
  ff.df = sp.Stream(ff.df);
  ff2.df = sp2.Stream(ff2.df);
  hwbb.UpdateNodes(f.df
    , g_is_modify_stream);
  hwbb2.UpdateNodes(f2.df
    , g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , g_is_modify_stream);
  hwbbsp2.UpdateNodes(ff2.df
    , g_is_modify_stream);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    auto left = n % g_nx == 0;
    auto right = n % g_nx == g_nx - 1;
    auto bottom = n / g_nx == 0;
    auto top = n / g_nx == g_ny - 1;
    // boundary on left & right
    if (bounce_back[n]) {
      CHECK_CLOSE(bottom || left ? f.df[n][0] : result_ne[f.df[n][0]],
          f.df[n][NE], zero_tol);
      CHECK_CLOSE(top || right ? f.df[n][0] : result_sw[f.df[n][0]],
          f.df[n][SW], zero_tol);
      CHECK_CLOSE(bottom || left ? ff.df[n][0] : result_ne[ff.df[n][0]],
          ff.df[n][NE], zero_tol);
      CHECK_CLOSE(top || right ? ff.df[n][0] : result_sw[ff.df[n][0]],
          ff.df[n][SW], zero_tol);
      CHECK_CLOSE(ff.df[n][0], ff.df[n][NW], zero_tol);
      CHECK_CLOSE(ff.df[n][0], ff.df[n][SE], zero_tol);
    }
    else {
      CHECK_CLOSE(bottom ? f.df[n][0] : result_ne[f.df[n][0]], f.df[n][NE],
          zero_tol);
      CHECK_CLOSE(top ? f.df[n][0] : result_sw[f.df[n][0]], f.df[n][SW],
          zero_tol);
      CHECK_CLOSE(result_ne[ff.df[n][0]], ff.df[n][NE], zero_tol);
      CHECK_CLOSE(result_sw[ff.df[n][0]], ff.df[n][SW], zero_tol);
      CHECK_CLOSE(bottom ? ff.df[n + height + 1][NW] : ff.df[n - g_nx + 1][NW],
          ff.df[n][NW], zero_tol);
      CHECK_CLOSE(top ? ff.df[n - height - 1][SE] : ff.df[n + g_nx - 1][SE],
          ff.df[n][SE], zero_tol);
    }
    CHECK_CLOSE(f.df[n][0], f.df[n][NW], zero_tol);
    CHECK_CLOSE(f.df[n][0], f.df[n][SE], zero_tol);
    // boundary on top & bottom
    if (bounce_back2[n]) {
      CHECK_CLOSE(bottom || left ? f2.df[n][0] : result_ne[f2.df[n][0]],
          f2.df[n][NE], zero_tol);
      CHECK_CLOSE(top || right ? f2.df[n][0] : result_sw[f2.df[n][0]],
          f2.df[n][SW], zero_tol);
      CHECK_CLOSE(bottom || left ? ff2.df[n][0] : result_ne[ff2.df[n][0]],
          ff2.df[n][NE], zero_tol);
      CHECK_CLOSE(top || right ? ff2.df[n][0] : result_sw[ff2.df[n][0]],
          ff2.df[n][SW], zero_tol);
    }
    else {
      CHECK_CLOSE(left ? f2.df[n][0] : result_ne[f2.df[n][0]], f2.df[n][NE],
          zero_tol);
      CHECK_CLOSE(right ? f2.df[n][0] : result_sw[f2.df[n][0]], f2.df[n][SW],
          zero_tol);
      CHECK_CLOSE(left ? ff2.df[n][0] : result_ne[ff2.df[n][0]], ff2.df[n][NE],
          zero_tol);
      CHECK_CLOSE(right ? ff2.df[n][0] : result_sw[ff2.df[n][0]], ff2.df[n][SW],
          zero_tol);
      CHECK_CLOSE(right ? ff2.df[n - width - g_nx][NW] : ff2.df[n][0],
          ff2.df[n][NW], zero_tol);
      CHECK_CLOSE(left ? ff2.df[n + width + g_nx][SE] : ff2.df[n][0],
          ff2.df[n][SE], zero_tol);
    }
    CHECK_CLOSE(f2.df[n][0], f2.df[n][NW], zero_tol);
    CHECK_CLOSE(f2.df[n][0], f2.df[n][SE], zero_tol);
  }  // n
}

TEST(StreamDiagonalNWSEOnGridBounceback)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  StreamD2Q9 sd(lm);
  StreamD2Q9 sd2(lm);
  StreamPeriodic sp(lm);
  StreamPeriodic sp2(lm);
  BouncebackNodes hwbb(lm
    , &sd);
  BouncebackNodes hwbb2(lm
    , &sd2);
  BouncebackNodes hwbbsp(lm
    , &sp);
  BouncebackNodes hwbbsp2(lm
    , &sp2);
  LatticeBoltzmann f(lm
    , ns
    , sd);
  LatticeBoltzmann f2(lm
    , ns
    , sd2);
  LatticeBoltzmann ff(lm
    , ns
    , sp);
  LatticeBoltzmann ff2(lm
    , ns
    , sp2);
  // same results for half-way bounceback with StreamD2Q9 on top & bottom and
  // left & right
  std::vector<bool> bounce_back(g_nx * g_ny, false);
  for (auto y = 0u; y < g_ny; ++y) {
    hwbb.AddNode(0, y);
    hwbb.AddNode(g_nx - 1, y);
    hwbbsp.AddNode(0, y);
    hwbbsp.AddNode(g_nx - 1, y);
    bounce_back[y * g_nx] = true;
    bounce_back[y * g_nx + g_nx - 1] = true;
  }
  std::vector<bool> bounce_back2(g_nx * g_ny, false);
  for (auto x = 0u; x < g_nx; ++x) {
    hwbb2.AddNode(x, 0);
    hwbb2.AddNode(x, g_ny - 1);
    hwbbsp2.AddNode(x, 0);
    hwbbsp2.AddNode(x, g_ny - 1);
    bounce_back2[x] = true;
    bounce_back2[(g_ny - 1) * g_nx + x] = true;
  }
  auto width = g_nx - 1;
  auto height = (g_ny - 1) * g_nx;
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
      f2.df[n] = nodes[(2 - x % 3 + (y + 2) % 3) % 3];
      ff.df[n] = nodes[(2 - x % 3 + (y + 2) % 3) % 3];
      ff2.df[n] = nodes[(2 - x % 3 + (y + 2) % 3) % 3];
    }  // x
  }  // y
  hwbb.UpdateNodes(f.df
    , !g_is_modify_stream);
  hwbb2.UpdateNodes(f2.df
    , !g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , !g_is_modify_stream);
  hwbbsp2.UpdateNodes(ff2.df
    , !g_is_modify_stream);
  f.df = sd.Stream(f.df);
  f2.df = sd2.Stream(f2.df);
  ff.df = sp.Stream(ff.df);
  ff2.df = sp2.Stream(ff2.df);
  hwbb.UpdateNodes(f.df
    , g_is_modify_stream);
  hwbb2.UpdateNodes(f2.df
    , g_is_modify_stream);
  hwbbsp.UpdateNodes(ff.df
    , g_is_modify_stream);
  hwbbsp2.UpdateNodes(ff2.df
    , g_is_modify_stream);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    auto left = n % g_nx == 0;
    auto right = n % g_nx == g_nx - 1;
    auto bottom = n / g_nx == 0;
    auto top = n / g_nx == g_ny - 1;

    // boundary on left & right
    if (bounce_back[n]) {
      CHECK_CLOSE(bottom || right ? f.df[n][0] : result_nw[f.df[n][0]],
          f.df[n][NW], zero_tol);
      CHECK_CLOSE(top || left ? f.df[n][0] : result_se[f.df[n][0]],
          f.df[n][SE], zero_tol);
      CHECK_CLOSE(bottom || right ? ff.df[n][0] : result_nw[ff.df[n][0]],
          ff.df[n][NW], zero_tol);
      CHECK_CLOSE(top || left ? ff.df[n][0] : result_se[ff.df[n][0]],
          ff.df[n][SE], zero_tol);
      CHECK_CLOSE(ff.df[n][0], ff.df[n][NE], zero_tol);
      CHECK_CLOSE(ff.df[n][0], ff.df[n][SW], zero_tol);
    }
    else {
      CHECK_CLOSE(bottom ? f.df[n][0] : result_nw[f.df[n][0]], f.df[n][NW],
          zero_tol);
      CHECK_CLOSE(top ? f.df[n][0] : result_se[f.df[n][0]], f.df[n][SE],
          zero_tol);
      CHECK_CLOSE(result_nw[ff.df[n][0]], ff.df[n][NW], zero_tol);
      CHECK_CLOSE(result_se[ff.df[n][0]], ff.df[n][SE], zero_tol);
      CHECK_CLOSE(bottom ? ff.df[n + height - 1][NE] : ff.df[n - g_nx - 1][NE],
          ff.df[n][NE], zero_tol);
      CHECK_CLOSE(top ? ff.df[n - height + 1][SW] : ff.df[n + g_nx + 1][SW],
          ff.df[n][SW], zero_tol);
    }
    CHECK_CLOSE(f.df[n][0], f.df[n][NE], zero_tol);
    CHECK_CLOSE(f.df[n][0], f.df[n][SW], zero_tol);
    // boundary on top & bottom
    if (bounce_back2[n]) {
      CHECK_CLOSE(bottom || right ? f2.df[n][0] : result_nw[f2.df[n][0]],
          f2.df[n][NW], zero_tol);
      CHECK_CLOSE(top || left ? f2.df[n][0] : result_se[f2.df[n][0]],
          f2.df[n][SE], zero_tol);
      CHECK_CLOSE(bottom || right ? ff2.df[n][0] : result_nw[ff2.df[n][0]],
          ff2.df[n][NW], zero_tol);
      CHECK_CLOSE(top || left ? ff2.df[n][0] : result_se[ff2.df[n][0]],
          ff2.df[n][SE], zero_tol);
    }
    else {
      CHECK_CLOSE(right ? f2.df[n][0] : result_nw[f2.df[n][0]], f2.df[n][NW],
          zero_tol);
      CHECK_CLOSE(left ? f2.df[n][0] : result_se[f2.df[n][0]], f2.df[n][SE],
          zero_tol);
      CHECK_CLOSE(right ? ff2.df[n][0] : result_nw[ff2.df[n][0]], ff2.df[n][NW],
          zero_tol);
      CHECK_CLOSE(left ? ff2.df[n][0] : result_se[ff2.df[n][0]], ff2.df[n][SE],
          zero_tol);
      CHECK_CLOSE(left ? ff2.df[n + width - g_nx][NE] : ff2.df[n][0],
          ff2.df[n][NE], zero_tol);
      CHECK_CLOSE(right ? ff2.df[n - width + g_nx][SW] : ff2.df[n][0],
          ff2.df[n][SW], zero_tol);
    }
    CHECK_CLOSE(f2.df[n][0], f2.df[n][NE], zero_tol);
    CHECK_CLOSE(f2.df[n][0], f2.df[n][SW], zero_tol);
  }  // n
}

TEST(InstantSourceToggle)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  StreamPeriodic sp(lm);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g
    , g_is_instant);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  cd.Collide(g.df);
  // check cd source is zero
  for (auto node : cd.source) CHECK_CLOSE(0.0, node, zero_tol);
}

TEST(InstantPointer)
{
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  StreamPeriodic sp(lm);
  CollisionCD cd(lm
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g
    , g_is_instant);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  g.TakeStep();
  // check cd source is zero
  for (auto node : cd.source) CHECK_CLOSE(0.0, node, zero_tol);
}

TEST(InterpolationStencils)
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
  ImmersedBoundaryMethod ibm(2
    , nsf.source
    , lm);
  for (auto x = -2.0; x <= 2.0; x += 0.01) {
    auto x_abs = fabs(x);
    CHECK_CLOSE(x_abs <= g_dx ? 1.0 - x_abs / g_dx : 0.0,
        ibm.Phi2(x, g_dx), loose_tol);
    if (x_abs <= 0.5) {
      CHECK_CLOSE((1.0 + sqrt(1.0 - 3.0 * x * x)) / 3.0, ibm.Phi3(x, g_dx),
          loose_tol);
    }
    else if (x_abs <= 1.5) {
      CHECK_CLOSE((5.0 - 3.0 * x_abs - sqrt(-2.0 + 6.0 * x_abs - 3.0 * x * x)) /
          6.0, ibm.Phi3(x, g_dx), loose_tol);
    }
    else {
      CHECK_CLOSE(0.0, ibm.Phi3(x, g_dx), loose_tol);
    }
    if (x_abs <= 1.0) {
      CHECK_CLOSE((3.0 - 2.0 * x_abs + sqrt(1.0 + 4.0 * x_abs - 4.0 * x * x)) /
          8.0, ibm.Phi4(x, g_dx), loose_tol);
    }
    else if (x_abs <= 2.0) {
      CHECK_CLOSE((5.0 - 2.0 * x_abs - sqrt(-7.0 + 12.0 * x_abs - 4.0 * x *
          x)) / 8.0, ibm.Phi4(x, g_dx), loose_tol);
    }
    else {
      CHECK_CLOSE(0.0, ibm.Phi4(x, g_dx), loose_tol);
    }
  }
}

TEST(CreateParticle)
{
  auto stiffness = -1.0;
  auto center_x = 10.0 * g_dx;
  auto center_y = 5.0 * g_dx;
  auto num_nodes = 10;
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  ParticleRigid cylinder(stiffness
    , num_nodes
    , center_x
    , center_y
    , lm);
  CHECK_CLOSE(center_x, cylinder.center.coord[0], zero_tol);
  CHECK_CLOSE(center_y, cylinder.center.coord[1], zero_tol);
  CHECK_CLOSE(center_x, cylinder.center.coord_ref[0], zero_tol);
  CHECK_CLOSE(center_y, cylinder.center.coord_ref[1], zero_tol);
  for (auto i = 0u; i < 2; ++i) {
    CHECK_CLOSE(0.0, cylinder.center.u[i], zero_tol);
    CHECK_CLOSE(0.0, cylinder.center.force[i], zero_tol);
  }  // i
  std::vector<double> x = {1.1 * g_dx, 2.1 * g_dx};
  std::vector<double> y = {1.2 * g_dx, 2.2 * g_dx};
  std::vector<double> x_ref = {1.3 * g_dx, 2.3 * g_dx};
  std::vector<double> y_ref = {1.4 * g_dx, 2.4 * g_dx};
  std::vector<double> u_x = {1.5, 2.5};
  std::vector<double> u_y = {1.6, 2.6};
  std::vector<double> force_x = {1.7, 2.7};
  std::vector<double> force_y = {1.8, 2.8};
  for (auto i = 0u; i < 2; ++i) {
    cylinder.AddNode(x[i]
    , y[i]
    , x_ref[i]
    , y_ref[i]
    , u_x[i]
    , u_y[i]
    , force_x[i]
    , force_y[i]);
  }  // i
  for (auto i = 0u; i < 2; ++i) {
    CHECK_CLOSE(x[i], cylinder.nodes[i].coord[0], zero_tol);
    CHECK_CLOSE(y[i], cylinder.nodes[i].coord[1], zero_tol);
    CHECK_CLOSE(x_ref[i], cylinder.nodes[i].coord_ref[0], zero_tol);
    CHECK_CLOSE(y_ref[i], cylinder.nodes[i].coord_ref[1], zero_tol);
    CHECK_CLOSE(u_x[i], cylinder.nodes[i].u[0], zero_tol);
    CHECK_CLOSE(u_y[i], cylinder.nodes[i].u[1], zero_tol);
    CHECK_CLOSE(force_x[i], cylinder.nodes[i].force[0], zero_tol);
    CHECK_CLOSE(force_y[i], cylinder.nodes[i].force[1], zero_tol);
  }  // i
}

TEST(CreateCylinderParticle)
{
  std::size_t num_nodes = 36;
  auto radius = 2.0;
  auto stiffness = -1.0;
  auto center_x = 11.0 * g_dx;
  auto center_y = 11.0 * g_dx;
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  ParticleRigid cylinder(stiffness
    , num_nodes
    , center_x
    , center_y
    , lm);
  cylinder.CreateCylinder(radius);
  for (auto node : cylinder.nodes) {
    auto x = node.coord[0] - center_x;
    auto y = node.coord[1] - center_y;
    // check is radius in lattice units
    CHECK_CLOSE(radius, sqrt(x * x + y * y), loose_tol);
  }  // node
}

TEST(ImmersedBoundaryForceReferencing)
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
  ImmersedBoundaryMethod ibm(2
    , nsf.source
    , lm);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    CHECK_CLOSE(nsf.source[n][0], ibm.fluid_force[n][0], zero_tol);
    CHECK_CLOSE(nsf.source[n][1], ibm.fluid_force[n][1], zero_tol);
  }  // n
}

TEST(ImmersedBoundaryClearVelocityForInterpolation)
{
  std::vector<double> u0 = {1.1, 1.2};
  std::size_t num_nodes = 36;
  auto radius = 2.0 * g_dx;  // 2 lattice units
  auto stiffness = -1.0;
  auto center_x = 11.0 * g_dx;
  auto center_y = 11.0 * g_dx;
  auto scaling = g_dx / g_dt;
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , u0);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  ParticleRigid cylinder(stiffness
    , num_nodes
    , center_x
    , center_y
    , lm);
  cylinder.CreateCylinder(radius);
  ImmersedBoundaryMethod ibm(2
    , nsf.source
    , lm);
  ibm.AddParticle(&cylinder);
  for (auto n = 0u; n < num_nodes; ++n) cylinder.nodes[n].u = {1.0, 1.1};
  ibm.InterpolateFluidVelocity();
  for (auto n = 0u; n < num_nodes; ++n) {
    CHECK_CLOSE(u0[0] / scaling, cylinder.nodes[n].u[0], loose_tol);
    CHECK_CLOSE(u0[1] / scaling, cylinder.nodes[n].u[1], loose_tol);
  }  // n
}

TEST(ImmersedBoundarySpreadForce)
{
  auto nx = 30u;
  auto ny = 20u;
  std::size_t num_nodes = 36;
  auto radius = 5.0 * g_dx;
  auto stiffness = -1.0;
  auto center_x = 11.0 * g_dx;
  auto center_y = 11.0 * g_dx;
  std::vector<std::size_t> cylinder_index = {188, 189, 190, 191, 192, 193, 194,
      217, 218, 219, 220, 221, 222, 223, 224, 225, 246, 247, 248, 254, 255, 256,
      276, 277, 285, 286, 306, 307, 315, 316, 336, 337, 345, 346, 347, 366, 367,
      375, 376, 377, 396, 397, 405, 406, 426, 427, 428, 434, 435, 436, 457, 458,
      459, 460, 461, 462, 463, 464, 465, 488, 489, 490, 491, 492, 493, 494, 521,
      522};
  std::vector<double> x_ans = {5737.607563, 22977.842969, 34925.881401,
      43224.096066, 34925.882030, 22977.842816, 5737.607608, 12347.872719,
      35589.530966, 19085.257920, 5331.070970, 695.796821, 5331.071022,
      19085.258165, 35589.530244, 12347.872719, 5737.607619, 35589.530063,
      9277.851160, 9277.851160, 35589.530786, 5737.607574, 22977.842778,
      19085.258226, 19085.257982, 22977.842930, 34925.882343, 5331.071035,
      5331.070983, 34925.881402, 43224.095754, 695.796821, 695.796821,
      43224.096378, 0.000000, 34925.881400, 5331.070957, 5331.071009,
      34925.881717, 0.000000, 22977.843007, 19085.257859, 19085.258104,
      22977.842854, 5737.607552, 35589.531147, 9277.851160, 9277.851160,
      35589.530424, 5737.607596, 12347.872719, 35589.529883, 19085.258287,
      5331.071048, 695.796757, 5331.070996, 19085.258043, 35589.530605,
      12347.872719, 5737.607630, 22977.842739, 34925.881408, 43224.095526,
      34925.881403, 22977.842892, 5737.607585, 0.000000, 0.000000};
  std::vector<double> y_ans = {6259.208250, 25066.737784, 38100.961528,
      47153.559345, 38100.962214, 25066.737617, 6259.208299, 13470.406603,
      38824.942872, 20820.281368, 5815.713786, 759.051078, 5815.713842,
      20820.281635, 38824.942084, 13470.406603, 6259.208311, 38824.941887,
      10121.292174, 10121.292174, 38824.942675, 6259.208262, 25066.737576,
      20820.281701, 20820.281434, 25066.737742, 38100.962556, 5815.713856,
      5815.713800, 38100.961530, 47153.559004, 759.051078, 759.051078,
      47153.559685, 0.000000, 38100.961527, 5815.713772, 5815.713828,
      38100.961873, 0.000000, 25066.737826, 20820.281301, 20820.281568,
      25066.737659, 6259.208238, 38824.943069, 10121.292174, 10121.292174,
      38824.942281, 6259.208287, 13470.406603, 38824.941690, 20820.281768,
      5815.713870, 759.051007, 5815.713814, 20820.281501, 38824.942478,
      13470.406603, 6259.208324, 25066.737534, 38100.961536, 47153.558755,
      38100.961531, 25066.737701, 6259.208275, 0.000000, 0.000000};
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  ParticleRigid cylinder(stiffness
    , num_nodes
    , center_x
    , center_y
    , lm);
  cylinder.CreateCylinder(radius);
  for (auto &node : cylinder.nodes) {
    node.force = {1.1, 1.2};
  }
  ImmersedBoundaryMethod ibm(2
    , nsf.source
    , lm);
  ibm.AddParticle(&cylinder);
  ibm.SpreadForce();
  auto i = 0;
  for (auto n = 0u; n < nx * ny; ++n) {
    if (n != cylinder_index[i]) {
      // checks that none IBM nodes are not affected by spread force
      CHECK_CLOSE(0.0, nsf.source[n][0], loose_tol);
      CHECK_CLOSE(0.0, nsf.source[n][1], loose_tol);
    }
    else {
      CHECK_CLOSE(x_ans[i], nsf.source[n][0], loose_tol);
      CHECK_CLOSE(y_ans[i], nsf.source[n][1], loose_tol);
      ++i;
    }
  }  // n
}

TEST(ImmersedBoundaryUpdateParticlePosition)
{
  auto nx = 30u;
  auto ny = 20u;
  std::size_t num_nodes = 36;
  auto radius = 2.0 * g_dx;
  auto stiffness = -1.0;
  auto center_x = 11.0 * g_dx;
  auto center_y = 11.0 * g_dx;
  auto particle_u = 1.0;
  auto particle_v = 1.1;
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  ParticleRigid cylinder(stiffness
    , num_nodes
    , center_x
    , center_y
    , lm);
  cylinder.CreateCylinder(radius);
  ImmersedBoundaryMethod ibm(2
    , nsf.source
    , lm);
  std::vector<std::vector<double>> exp_coord;
  for (auto &node : cylinder.nodes) {
    exp_coord.push_back({node.coord[0] + particle_u * g_dt, node.coord[1] +
        particle_v * g_dt});
    node.u = {particle_u, particle_v};
  }
  ibm.AddParticle(&cylinder);
  ibm.UpdateParticlePosition();
  for (auto n = 0u; n < num_nodes; ++n) {
    CHECK_CLOSE(exp_coord[n][0], cylinder.nodes[n].coord[0], loose_tol);
    CHECK_CLOSE(exp_coord[n][1], cylinder.nodes[n].coord[1], loose_tol);
  }  // n
}

TEST(RigidParticleComputeParticleForce)
{
  std::size_t num_nodes = 36;
  auto radius = 2.0 * g_dx;
  auto stiffness = -1.0;
  auto center_x = 11.0 * g_dx;
  auto center_y = 11.0 * g_dx;
  auto displacement = 0.1;
  LatticeD2Q9 lm(g_ny
    , g_nx
    , g_dx
    , g_dt
    , g_u0);
  ParticleRigid cylinder(stiffness
    , num_nodes
    , center_x
    , center_y
    , lm);
  cylinder.CreateCylinder(radius);
  auto force = -stiffness * 2.0 * g_pi * radius / num_nodes * displacement;
  for (auto &node : cylinder.nodes) {
    node.coord[0] += displacement;
    node.coord[1] += displacement;
  }  // node
  cylinder.ComputeForces();
  for (auto node : cylinder.nodes) {
    CHECK_CLOSE(force, node.force[0], loose_tol);
    CHECK_CLOSE(force, node.force[1], loose_tol);
  }  // node
}

TEST(MobileRigidParticle)
{
  auto nx = 30u;
  auto ny = 20u;
  std::size_t num_nodes = 36;
  auto radius = 2.0 * g_dx;
  auto stiffness = -1.0;
  auto center_x = 11.0 * g_dx;
  auto center_y = 11.0 * g_dx;
  auto particle_u = 1.0;
  auto particle_v = 1.1;
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , g_u0);
  CollisionNSF nsf(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f);
  ParticleRigid cylinder(stiffness
    , num_nodes
    , center_x
    , center_y
    , lm);
  cylinder.CreateCylinder(radius);
  cylinder.ChangeMobility(true);
  ImmersedBoundaryMethod ibm(2
    , nsf.source
    , lm);
  std::vector<std::vector<double>> exp_coord;
  for (auto &node : cylinder.nodes) {
    exp_coord.push_back({node.coord[0] + particle_u * g_dt, node.coord[1] +
        particle_v * g_dt});
    node.u = {particle_u, particle_v};
  }
  ibm.AddParticle(&cylinder);
  ibm.UpdateParticlePosition();
  for (auto n = 0u; n < num_nodes; ++n) {
    CHECK_CLOSE(exp_coord[n][0], cylinder.nodes[n].coord[0], loose_tol);
    CHECK_CLOSE(exp_coord[n][1], cylinder.nodes[n].coord[1], loose_tol);
    CHECK_CLOSE(exp_coord[n][0], cylinder.nodes[n].coord_ref[0], loose_tol);
    CHECK_CLOSE(exp_coord[n][1], cylinder.nodes[n].coord_ref[1], loose_tol);
  }  // n
}
}

SUITE(TestSteadyState)
{
static const double zero_tol = 1e-20;
static const double loose_tol = 1e-5;
static const double g_dx = 0.0316;
static const double g_dt = 0.001;
static const double g_d_coeff = 0.2;
static const double g_k_visco = 0.2;
static const double g_rho0_f = 1.1;
static const double g_rho0_g = 1.2;

TEST(TerminationConditionSteadyState)
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
  auto u_prev = lm.u;
  auto t = 0;
  // this simulation is known to reach steady state in 1000 < t < 1500 ts
  for (auto t = 0; t < 1000; ++t) {
    u_prev = lm.u;
    f.TakeStep();
    CHECK_EQUAL(false, CheckSteadyState(u_prev, lm.u, loose_tol));
  }
  for (auto t = 0; t < 500; ++t) {
    u_prev = lm.u;
    f.TakeStep();
  }
  CHECK_EQUAL(true, CheckSteadyState(u_prev, lm.u, loose_tol));
}

TEST(TerminationConditionSteadyStateZHInlet)
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
  std::vector<std::vector<double>> u_curr;
  auto u_prev = u_curr;
  for (auto n = nx; n < nx * (ny - 1); ++n) u_curr.push_back(lm.u[n]);
  auto t = 0;
  // this simulation is known to reach steady state in 1500 < t < 2000 ts
  for (auto t = 0; t < 1500; ++t) {
    u_prev = u_curr;
    f.TakeStep();
    auto i = 0;
    for (auto n = nx; n < nx * (ny - 1); ++n) u_curr[i++] = lm.u[n];
    CHECK_EQUAL(false, CheckSteadyState(u_prev, u_curr, loose_tol));
  }
  for (auto t = 0; t < 500; ++t) {
    u_prev = u_curr;
    f.TakeStep();
    auto i = 0;
    for (auto n = nx; n < nx * (ny - 1); ++n) u_curr[i++] = lm.u[n];
  }
  CHECK_EQUAL(true, CheckSteadyState(u_prev, u_curr, loose_tol));
}
}
