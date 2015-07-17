#include "CollisionModel.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Algorithm.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"
#include "Printing.hpp"

CollisionModel::CollisionModel(LatticeModel &lm)
  : f_eq {},
    g_eq {},
    rho_f {},
    rho_g {},
    is_ns {},
    is_cd {},
    lm_ (lm),
    c_ {lm.GetLatticeSpeed()}
{}

CollisionModel::CollisionModel(LatticeModel &lm
  , const std::vector<double> &initial_density)
  : f_eq {},
    g_eq {},
    rho_f {initial_density},
    rho_g {},
    is_ns {},
    is_cd {},
    lm_ (lm),
    c_ {lm.GetLatticeSpeed()}
{}

void CollisionModel::ComputeEq(std::vector<std::vector<double>> &lattice_eq
  , const std::vector<double> &rho)
{
  auto nc = lm_.GetNumberOfDirections();
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (auto n = 0u; n < nx * ny; ++n) {
    double u_sqr = InnerProduct(lm_.u[n], lm_.u[n]);
    u_sqr /= 2.0 * cs_sqr_;
    for (auto i = 0u; i < nc; ++i) {
      double c_dot_u = InnerProduct(lm_.e[i], lm_.u[n]);
      c_dot_u /= cs_sqr_;
      lattice_eq[n][i] = lm_.omega[i] * rho[n] * (1.0 + c_dot_u *
          (1.0 + c_dot_u / 2.0) - u_sqr);
    }  // i
  }  // n
}
