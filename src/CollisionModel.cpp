#include "CollisionModel.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

CollisionModel::CollisionModel(LatticeModel &lm
  , double initial_density)
  : edf {},
    rho {},
    lm_ (lm),
    tau_ {0},
    c_ {lm.GetLatticeSpeed()}
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nc = lm_.GetNumberOfDirections();
  auto lat_size = nx * ny;
  rho.assign(lat_size, initial_density);
  edf.assign(lat_size, std::vector<double>(nc, 0.0));
  ComputeEq();
}

void CollisionModel::ComputeEq()
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
      edf[n][i] = lm_.omega[i] * rho[n] * (1.0 + c_dot_u * (1.0 + c_dot_u /
          2.0) - u_sqr);
    }  // i
  }  // n
}
