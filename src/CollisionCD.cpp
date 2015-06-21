#include "CollisionCD.hpp"
#include <iostream>
#include <stdexcept>
#include "LatticeModel.hpp"

// specifies the base class constructor to call
// https://stackoverflow.com/questions/10282787/calling-the-base-class-
// constructor-in-the-derived-class-constructor
CollisionCD::CollisionCD(LatticeModel &lm
  , std::vector<std::vector<double>> &lat
  , const std::vector<std::vector<std::size_t>> &position
  , const std::vector<double> &strength
  , double diffusion_coefficient)
  : Collision(lm)
{
  lattice_eq_ = lat;
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
  source_.assign(nx * ny, 0.0);
  auto it_strength = begin(strength);
  for (auto pos : position) {
    if (pos.size() != nd)
        throw std::runtime_error("Position doesn't match dimensions");
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    source_[pos[1] * nx + pos[0]] = *it_strength++;
  }  // pos
}

void CollisionCD::ApplyForce()
{
//  std::cout << source_[0] <<std::endl;
  std::cout << lm_.GetNumberOfColumns() << std::endl;
}
