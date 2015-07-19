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
  auto dt = lm.GetTimeStep();
  // tau_ formula from "Discrete lattice effects on the forcing term in
  // the lattice Boltzmann method" Guo2002
  tau_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
}

void CollisionNS::Collide(std::vector<std::vector<double>> &lattice)
{
  auto nc = lm_.GetNumberOfDirections();
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (auto n = 0u; n < nx * ny; ++n) {
    for (auto i = 0u; i < nc; ++i) {
      lattice[n][i] += (edf[n][i] - lattice[n][i]) / tau_;
    }  // i
  }  // n
}
