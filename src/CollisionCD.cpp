#include "CollisionCD.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

// specifies the base class constructor to call
// https://stackoverflow.com/questions/10282787/calling-the-base-class-
// constructor-in-the-derived-class-constructor
CollisionCD::CollisionCD(LatticeModel &lm
  , const std::vector<std::vector<std::size_t>> &source_position
  , const std::vector<double> &source_strength
  , double diffusion_coefficient
  , double initial_density_g
  , bool is_instant)
  : CollisionModel(lm, initial_density_g),
    source {},
    is_instant_ {is_instant}
{
  const auto dt = lm.GetTimeStep();
  // tau_ formula from "A new scheme for source term in LBGK model for
  // convection diffusion equation"
  tau_ = 0.5 + diffusion_coefficient / cs_sqr_ / dt;
  CollisionCD::InitSource(source_position, source_strength);
}

void CollisionCD::InitSource(
    const std::vector<std::vector<std::size_t>> &source_position
  , const std::vector<double> &source_strength)
{
  if (source_position.size() != source_strength.size())
      throw std::runtime_error("Insufficient source information");
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  const auto nd = lm_.GetNumberOfDimensions();
  source.assign(nx * ny, 0.0);
  auto it_strength = begin(source_strength);
  for (auto pos : source_position) {
    if (pos.size() != nd) throw std::runtime_error("Dimensions mismatch");
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    source[pos[1] * nx + pos[0]] = *it_strength++;
  }  // pos
}

void CollisionCD::ComputeMacroscopicProperties(
      const std::vector<std::vector<double>> &df)
{
  rho = CollisionCD::ComputeRho(df);
}

void CollisionCD::Collide(std::vector<std::vector<double>> &df)
{
  const auto nc = lm_.GetNumberOfDirections();
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  const auto dt = lm_.GetTimeStep();
  for (auto n = 0u; n < nx * ny; ++n) {
    if (!skip[n]) {
      for (auto i = 0u; i < nc; ++i) {
        double c_dot_u = InnerProduct(lm_.e[i], lm_.u[n]);
        c_dot_u /= cs_sqr_;
        // source term using forward scheme, theta = 0
        const auto src_i = lm_.omega[i] * source[n] * (1.0 + (1.0 - 0.5 /
            tau_) * c_dot_u);
        df[n][i] += (edf[n][i] - df[n][i]) / tau_ + dt * src_i;
      }  // i
    }
  }  // n
  if (is_instant_) CollisionCD::KillSource();
}

void CollisionCD::KillSource()
{
//  for (auto &node : source) node = 0.0;
  source = std::vector<double>(source.size(), 0.0);
}
