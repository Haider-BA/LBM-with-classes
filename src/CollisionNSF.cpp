#include "CollisionNSF.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

// specifies the base class constructor to call
// https://stackoverflow.com/questions/10282787/calling-the-base-class-
// constructor-in-the-derived-class-constructor
CollisionNSF::CollisionNSF(LatticeModel &lm
  , const std::vector<std::vector<std::size_t>> &position
  , const std::vector<std::vector<double>> &strength
  , double kinematic_viscosity
  , double initial_density_f)
  : CollisionNS(lm, kinematic_viscosity, initial_density_f),
    source {}
{
  // tau_ already initialized with calling of CollisionNS
//  InitSource(position, strength);
}

void CollisionNSF::CollideNS(std::vector<std::vector<double>> &lattice)
{
  auto nc = lm_.GetNumberOfDirections();
  auto nd = lm_.GetNumberOfDimensions();
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto dt = lm_.GetTimeStep();
  for (auto y = 0u; y < ny; ++y) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      for (auto i = 0u; i < nc; ++i) {
        double c_dot_u = InnerProduct(lm_.e[i], lm_.u[n]);
        c_dot_u /= cs_sqr_;
        double src_dot_product = 0.0;
        for (auto d = 0u; d < nd; ++d) {
          src_dot_product += (lm_.e[i][d] - lm_.u[n][d] + c_dot_u *
              lm_.e[i][d]) * source[n][d];
        }  // d
        src_dot_product /= cs_sqr_ / rho_f[n];
        auto src_i = (1.0 - 0.5 / tau_f_) * lm_.omega[i] * src_dot_product;
        lattice[n][i] += (f_eq[n][i] - lattice[n][i]) / tau_f_ +
            dt * src_i;
      }  // i
    }  // x
  }  // y
}
