#include "CollisionCD.hpp"
#include <iostream>
#include <stdexcept>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

// specifies the base class constructor to call
// https://stackoverflow.com/questions/10282787/calling-the-base-class-
// constructor-in-the-derived-class-constructor
CollisionCD::CollisionCD(LatticeModel &lm
  , const std::vector<std::vector<std::size_t>> &position
  , const std::vector<double> &strength
  , double diffusion_coefficient
  , double initial_density_g)
  : Collision(lm, initial_density_g),
    source {}
{
  auto dt = lm.GetTimeStep();
  // tau_ formula from "A new scheme for source term in LBGK model for
  // convection diffusion equation"
  tau_ = 0.5 + diffusion_coefficient / cs_sqr_ / dt;
  InitSource(position, strength);
}

void CollisionCD::InitSource(
    const std::vector<std::vector<std::size_t>> &position
  , const std::vector<double> &strength)
{
  if (position.size() != strength.size())
      throw std::runtime_error("Insufficient source information");
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nd = lm_.GetNumberOfDimensions();
  source.assign(nx * ny, 0.0);
  auto it_strength = begin(strength);
  for (auto pos : position) {
    if (pos.size() != nd)
        throw std::runtime_error("Position doesn't match dimensions");
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    source[pos[1] * nx + pos[0]] = *it_strength++;
  }  // pos
}


void CollisionCD::Collide(std::vector<std::vector<double>> &lattice)
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
      auto src_i = lm_.omega[i] * source[n] * (1.0 + (1.0 - 0.5 / tau_) *
          c_dot_u);
      lattice[n][i] += (lattice_eq[n][i] - lattice[n][i]) / tau_ + dt * src_i;
    }  // i
  }  // n
}

void CollisionCD::KillSource()
{
  for (auto &node : source) node = 0.0;
}
