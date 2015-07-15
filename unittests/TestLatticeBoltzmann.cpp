#include <iostream>
#include <stdexcept>  // runtime_error
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"
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

TEST(test)
{
  LatticeD2Q9 lm(g_ny, g_nx, g_dx, g_dt, g_u0);
  CollisionCD cd;
  CollisionNS ns;
  LatticeBoltzmann lbm(0.1, lm, ns, cd);
}
}
