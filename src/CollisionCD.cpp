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
  , double initial_density_g)
  : CollisionModel(lm),
    source {},
    tau_g_ {0.0}
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nc = lm_.GetNumberOfDirections();
  auto dt = lm.GetTimeStep();
  auto lattice_size = nx * ny;
  rho_g.assign(lattice_size, initial_density_g);
  g_eq.assign(lattice_size, std::vector<double>(nc, 0.0));
  CollisionCD::ComputeEq(g_eq, rho_g);
  // tau_ formula from "A new scheme for source term in LBGK model for
  // convection diffusion equation"
  tau_g_ = 0.5 + diffusion_coefficient / cs_sqr_ / dt;
  InitSource(source_position, source_strength);
}

void CollisionCD::InitSource(
    const std::vector<std::vector<std::size_t>> &source_position
  , const std::vector<double> &source_strength)
{
  if (source_position.size() != source_strength.size())
      throw std::runtime_error("Insufficient source information");
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nd = lm_.GetNumberOfDimensions();
  source.assign(nx * ny, 0.0);
  auto it_strength = begin(source_strength);
  for (auto pos : source_position) {
    if (pos.size() != nd)
        throw std::runtime_error("Position doesn't match dimensions");
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    source[pos[1] * nx + pos[0]] = *it_strength++;
  }  // pos
}


void CollisionCD::CollideCD(std::vector<std::vector<double>> &lattice)
{
  auto nc = lm_.GetNumberOfDirections();
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto dt = lm_.GetTimeStep();
  for (auto n = 0u; n < nx * ny; ++n) {
    for (auto i = 0u; i < nc; ++i) {
      double c_dot_u = InnerProduct(lm_.e[i], lm_.u[n]);
      c_dot_u /= cs_sqr_;
      // source term using forward scheme, theta = 0
      auto src_i = lm_.omega[i] * source[n] * (1.0 + (1.0 - 0.5 / tau_g_) *
          c_dot_u);
      lattice[n][i] += (g_eq[n][i] - lattice[n][i]) / tau_g_ + dt * src_i;
    }  // i
  }  // n
}

void CollisionCD::KillSource()
{
  for (auto &node : source) node = 0.0;
}
