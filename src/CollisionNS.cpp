#include "CollisionNS.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"
#include "Printing.hpp"

// specifies the base class constructor to call
// https://stackoverflow.com/questions/10282787/calling-the-base-class-
// constructor-in-the-derived-class-constructor
CollisionNS::CollisionNS(LatticeModel &lm
  , double kinematic_viscosity
  , double initial_density_f)
  : CollisionModel(lm),
    tau_f_ {0.0}
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nc = lm_.GetNumberOfDirections();
  auto dt = lm.GetTimeStep();
  auto lattice_size = nx * ny;
  rho_f.assign(lattice_size, initial_density_f);
  f_eq.assign(lattice_size, std::vector<double>(nc, 0.0));
  CollisionNS::ComputeEq(f_eq, rho_f);
  // tau_ formula from "Discrete lattice effects on the forcing term in
  // the lattice Boltzmann method" Guo2002
  tau_f_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
//  InitSource(position, strength);
  is_ns = true;
}

CollisionNS::CollisionNS(LatticeModel &lm
  , double kinematic_viscosity
  , const std::vector<double> &initial_density_f)
  : CollisionModel(lm),
    tau_f_ {0.0}
{
  auto dt = lm.GetTimeStep();
  // tau_ formula from "Discrete lattice effects on the forcing term in
  // the lattice Boltzmann method" Guo2002
//  tau_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
  tau_f_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
//  InitSource(position, strength);
}

//void CollisionNS::InitSource(
//    const std::vector<std::vector<std::size_t>> &position
//  , const std::vector<std::vector<double>> &strength)
//{
//  if (position.size() != strength.size())
//      throw std::runtime_error("Insufficient source information");
//  auto nx = lm_.GetNumberOfColumns();
//  auto ny = lm_.GetNumberOfRows();
//  auto nd = lm_.GetNumberOfDimensions();
//  source.assign(nx * ny, {0.0, 0.0});
//  auto it_strength = begin(strength);
//  for (auto pos : position) {
//    if (pos.size() != nd)
//        throw std::runtime_error("Position doesn't match dimensions");
//    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
//    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
//    source[pos[1] * nx + pos[0]] = *it_strength++;
//  }  // pos
//}

void CollisionNS::CollideNS(std::vector<std::vector<double>> &lattice)
{
  auto nc = lm_.GetNumberOfDirections();
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (auto n = 0u; n < nx * ny; ++n) {
    for (auto i = 0u; i < nc; ++i) {
      lattice[n][i] += (f_eq[n][i] - lattice[n][i]) / tau_f_;
    }  // i
  }  // n
}

void CollisionNS::CollideLid(std::vector<std::vector<double>> &lattice)
{
  auto nc = lm_.GetNumberOfDirections();
//  auto nd = lm_.GetNumberOfDimensions();
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
//  auto dt = lm_.GetTimeStep();
  for (auto y = 0u; y < ny; ++y) {
    for (auto x = 0u; x < nx; ++x) {
      if (y == ny - 1 && !x && x != nx - 1) continue;
      auto n = y * nx + x;
      for (auto i = 0u; i < nc; ++i) {
//        double c_dot_u = InnerProduct(lm_.e[i], lm_.u[n]);
//        c_dot_u /= cs_sqr_;
//        double src_dot_product = 0.0;
//        for (auto d = 0u; d < nd; ++d) {
//          src_dot_product += (lm_.e[i][d] - lm_.u[n][d] + c_dot_u *
//              lm_.e[i][d]) * source[n][d];
//        }  // d
//        src_dot_product /= cs_sqr_ / rho[n];
//        auto src_i = (1.0 - 0.5 / tau_) * lm_.omega[i] * src_dot_product;
        lattice[n][i] += (f_eq[n][i] - lattice[n][i]) / tau_f_;// + dt * src_i;
      }  // i
    }  // x
  }  // y
}

//void CollisionNS::KillSource()
//{
//  auto nd = lm_.GetNumberOfDimensions();
//  std::vector<double> zeros(nd, 0.0);
//  for (auto &node : source) node = zeros;
//}
