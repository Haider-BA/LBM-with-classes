#include <iostream>
#include <stdexcept>  // runtime_error
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "CollisionNSCD.hpp"
#include "CollisionNSF.hpp"
#include "CollisionNSFCD.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"
#include "Printing.hpp"
#include "UnitTest++.h"

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
static const std::size_t g_nx = 7;
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
  CollisionNSCD nscd(lm
    , g_k_visco
    , g_rho0_f
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  CollisionNSFCD nsfcd(lm
    , g_src_pos_f
    , g_src_str_f
    , g_k_visco
    , g_rho0_f
    , g_src_pos_g
    , g_src_str_g
    , g_d_coeff
    , g_rho0_g);
  for (auto n = 0u; n < g_nx * g_ny; ++n) {
    CHECK_CLOSE(g_rho0_f, ns.rho_f[n], zero_tol);
    CHECK_CLOSE(g_rho0_f, nsf.rho_f[n], zero_tol);
    CHECK_CLOSE(g_rho0_g, cd.rho_g[n], zero_tol);
    CHECK_CLOSE(g_rho0_f, nscd.rho_f[n], zero_tol);
    CHECK_CLOSE(g_rho0_g, nscd.rho_g[n], zero_tol);
    CHECK_CLOSE(g_rho0_f, nsfcd.rho_f[n], zero_tol);
    CHECK_CLOSE(g_rho0_g, nsfcd.rho_g[n], zero_tol);
  }
}
}
