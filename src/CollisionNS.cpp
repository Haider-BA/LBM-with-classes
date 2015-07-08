#include "CollisionNS.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

// specifies the base class constructor to call
// https://stackoverflow.com/questions/10282787/calling-the-base-class-
// constructor-in-the-derived-class-constructor
CollisionNS::CollisionNS(LatticeModel &lm
  , const std::vector<std::vector<std::size_t>> &position
  , const std::vector<std::vector<double>> &strength
  , double kinematic_viscosity
  , double initial_density_f)
  : Collision(lm, initial_density_f),
    source {}
{
  auto dt = lm.GetTimeStep();
  // tau_ formula from "Discrete lattice effects on the forcing term in
  // the lattice Boltzmann method" Guo2002
  tau_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
  InitSource(position, strength);
}

CollisionNS::CollisionNS(LatticeModel &lm
  , const std::vector<std::vector<std::size_t>> &position
  , const std::vector<std::vector<double>> &strength
  , double kinematic_viscosity
  , const std::vector<double> &initial_density_f)
  : Collision(lm, initial_density_f),
    source {}
{
  auto dt = lm.GetTimeStep();
  // tau_ formula from "Discrete lattice effects on the forcing term in
  // the lattice Boltzmann method" Guo2002
  tau_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
  InitSource(position, strength);
}

void CollisionNS::InitSource(
    const std::vector<std::vector<std::size_t>> &position
  , const std::vector<std::vector<double>> &strength)
{
  if (position.size() != strength.size())
      throw std::runtime_error("Insufficient source information");
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nd = lm_.GetNumberOfDimensions();
  if (source.empty()) source.assign(nx * ny, {0.0, 0.0});
  auto it_strength = begin(strength);
  for (auto pos : position) {
    if (pos.size() != nd)
        throw std::runtime_error("Position doesn't match dimensions");
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    auto temp = source[pos[1] * nx + pos[0]][0];
    source[pos[1] * nx + pos[0]] = *it_strength++;
  }  // pos
}

void CollisionNS::Collide(std::vector<std::vector<double>> &lattice)
{
  auto nc = lm_.GetNumberOfDirections();
  auto nd = lm_.GetNumberOfDimensions();
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto dt = lm_.GetTimeStep();
  for (auto n = 0u; n < nx * ny; ++n) {
    for (auto i = 0u; i < nc; ++i) {
      double c_dot_u = InnerProduct(lm_.e[i], lm_.u[n]);
      c_dot_u /= cs_sqr_;
      double src_dot_product = 0.0;
      for (auto d = 0u; d < nd; ++d) {
        src_dot_product += (lm_.e[i][d] - lm_.u[n][d] + c_dot_u * lm_.e[i][d]) *
            source[n][d];
      }  // d
      src_dot_product /= cs_sqr_ / rho[n];
      auto src_i = (1.0 - 0.5 / tau_) * lm_.omega[i] * src_dot_product;
      lattice[n][i] += (lattice_eq[n][i] - lattice[n][i]) / tau_ + dt * src_i;
    }  // i
  }  // n
}

void CollisionNS::KillSource()
{
  auto nd = lm_.GetNumberOfDimensions();
  std::vector<double> zeros(nd, 0.0);
  for (auto &node : source) node = zeros;
}
