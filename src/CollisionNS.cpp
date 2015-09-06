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
  , double kinematic_viscosity
  , double initial_density_f)
  : CollisionModel(lm, initial_density_f)
{
  const auto dt = lm.GetTimeStep();
  // tau_ formula from "Discrete lattice effects on the forcing term in
  // the lattice Boltzmann method" Guo2002
  tau_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
}

CollisionNS::CollisionNS(LatticeModel &lm
  , double kinematic_viscosity
  , const std::vector<double> &initial_density_f)
  : CollisionModel(lm, initial_density_f)
{
  const auto dt = lm.GetTimeStep();
  // tau_ formula from "Discrete lattice effects on the forcing term in
  // the lattice Boltzmann method" Guo2002
  tau_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
}

std::vector<std::vector<double>> CollisionNS::ComputeU(
    const std::vector<std::vector<double>> &df)
{
  std::vector<std::vector<double>> result;
  auto index = 0u;
  for (auto node : df) result.push_back(GetFirstMoment(node, lm_.e));
  for (auto &node : result) {
    for (auto &d : node) d /= rho[index];
    ++ index;
  }  // node
  return result;
}

void CollisionNS::ComputeMacroscopicProperties(
      const std::vector<std::vector<double>> &df)
{
  rho = CollisionNS::ComputeRho(df);
  lm_.u = CollisionNS::ComputeU(df);
}

void CollisionNS::Collide(std::vector<std::vector<double>> &lattice)
{
  const auto nc = lm_.GetNumberOfDirections();
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  for (auto n = 0u; n < nx * ny; ++n) {
    if (!skip[n]) {
      for (auto i = 0u; i < nc; ++i) {
        lattice[n][i] += (edf[n][i] - lattice[n][i]) / tau_;
      }  // i
    }
  }  // n
}
