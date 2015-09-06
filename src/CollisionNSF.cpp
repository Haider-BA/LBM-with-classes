#include "CollisionNSF.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

CollisionNSF::CollisionNSF(LatticeModel &lm
  , const std::vector<std::vector<std::size_t>> &source_position
  , const std::vector<std::vector<double>> &source_strength
  , double kinematic_viscosity
  , double initial_density_f)
  : CollisionNS(lm, kinematic_viscosity, initial_density_f),
    source {}
{
  const auto dt = lm.GetTimeStep();
  // tau_ formula from "Discrete lattice effects on the forcing term in
  // the lattice Boltzmann method" Guo2002
  tau_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
  CollisionNSF::InitSource(source_position, source_strength);
}

void CollisionNSF::InitSource(
    const std::vector<std::vector<std::size_t>> &source_position
  , const std::vector<std::vector<double>> &source_strength)
{
  if (source_position.size() != source_strength.size())
      throw std::runtime_error("Insufficient source information");
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  const auto nd = lm_.GetNumberOfDimensions();
  source.assign(nx * ny, std::vector<double>(nd, 0.0));
  auto it_strength = begin(source_strength);
  for (auto pos : source_position) {
    if (pos.size() != nd) throw std::runtime_error("Dimensions mismatch");
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    source[pos[1] * nx + pos[0]] = *it_strength++;
  }  // pos
}

std::vector<std::vector<double>> CollisionNSF::ComputeU(
    const std::vector<std::vector<double>> &df)
{
  const auto nd = lm_.GetNumberOfDimensions();
  const auto dt = lm_.GetTimeStep();
  std::vector<std::vector<double>> result;
  auto index = 0u;
  for (auto node : df) {
    result.push_back(GetFirstMoment(node, lm_.e));
    for (auto d = 0u; d < nd; ++d) {
      result[index][d] += 0.5 * dt * source[index][d] * rho[index];
      result[index][d] /= rho[index];
    }  // d
    ++index;
  }
  return result;
}

void CollisionNSF::ComputeMacroscopicProperties(
      const std::vector<std::vector<double>> &df)
{
  rho = CollisionNSF::ComputeRho(df);
  lm_.u = CollisionNSF::ComputeU(df);
}

void CollisionNSF::Collide(std::vector<std::vector<double>> &lattice)
{
  const auto nc = lm_.GetNumberOfDirections();
  const auto nd = lm_.GetNumberOfDimensions();
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  const auto dt = lm_.GetTimeStep();
  for (auto n = 0u; n < nx * ny; ++n) {
    if (!skip[n]) {
      for (auto i = 0u; i < nc; ++i) {
        double c_dot_u = InnerProduct(lm_.e[i], lm_.u[n]);
        c_dot_u /= cs_sqr_;
        double src_dot_product = 0.0;
        for (auto d = 0u; d < nd; ++d) {
          src_dot_product += (lm_.e[i][d] - lm_.u[n][d] + c_dot_u *
              lm_.e[i][d]) * source[n][d];
        }  // d
        src_dot_product /= cs_sqr_ / rho[n];
        auto src_i = (1.0 - 0.5 / tau_) * lm_.omega[i] * src_dot_product;
        lattice[n][i] += (edf[n][i] - lattice[n][i]) / tau_ + dt * src_i;
      }  // i
    }
  }  // n
}
