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
TEST(Constructor)
{
  std::size_t ny = 4;
  std::size_t nx = 5;
  double dx = 0.0316;
  double dt = 0.001;
  double t_total = 1.0;
  double d_coeff = 0.2;
  double k_visco = 0.2;
  double rho0_f = 1.0;
  double rho0_g = 1.0;
  std::vector<double> u0{1, 3};
  std::vector<std::vector<std::size_t>> src_pos_f = {{1, 1}};
  std::vector<std::vector<double>> src_str_f = {{1.5, 1.5}};
  std::vector<std::vector<std::size_t>> src_pos_g = {{2, 2}};
  std::vector<double> src_str_g = {2.5};
  std::vector<std::vector<std::size_t>> obs_pos;
  bool is_ns = true;
  bool is_cd = true;
  bool is_instant = true;
  bool no_obstacles = false;
  LatticeD2Q9 lm(ny, nx, dx, dt);
  LatticeBoltzmann lbm(t_total, d_coeff, k_visco,
      rho0_f, rho0_g, u0, src_pos_f, src_str_f, src_pos_g, src_str_g, obs_pos,
      is_ns, is_cd, is_instant, no_obstacles, lm);

  std::vector<std::vector<double>> lat = {{10, 11}, {12, 13}};
  CollisionCD ccd(lm, lat, src_pos_g, src_str_g, d_coeff);
  CollisionNS cns(lm, lat, src_pos_f, src_str_f, k_visco);
  ccd.Collide();
  ccd.ApplyForce();
  lbm.Print(ccd.source_);
  lbm.Print(1, cns.source_);
}
}
