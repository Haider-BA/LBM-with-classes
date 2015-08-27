#include "ImmersedBoundaryMethod.hpp"
#include <iostream>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

ImmersedBoundaryMethod::ImmersedBoundaryMethod(int interpolation_stencil
  , std::vector<std::vector<double>> &lattice_force
  , LatticeModel &lm)
  : particles{},
    fluid_force {lattice_force},
    interpolation_stencil_ {interpolation_stencil},
    lm_ {lm}
{}

void ImmersedBoundaryMethod::AddParticle(Particle* particle)
{
  particles.push_back(particle);
}

void ImmersedBoundaryMethod::InterpolateFluidVelocity()
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  // pointer always editable, do not need &
  for (auto particle : particles) {
    for (auto &node : particle -> nodes) {
      node.u = {0.0, 0.0};
      auto x_particle = node.coord[0];
      auto y_particle = node.coord[1];
      auto x_fluid = static_cast<std::size_t>(x_particle);
      auto y_fluid = static_cast<std::size_t>(y_particle);
      for (auto x = x_fluid; x < x_fluid + 2; ++x) {
        for (auto y = y_fluid; y < y_fluid + 2; ++y) {
          auto n = (y % ny) * nx + x % nx;
          auto weight = Dirac2(x_particle - x, y_particle - y);
          node.u[0] += lm_.u[n][0] * weight;
          node.u[1] += lm_.u[n][1] * weight;
        }  // y
      }  // x
    }  // node
  }  // particle
}

void ImmersedBoundaryMethod::SpreadForce()
{
  // The two-point interpolation stencil (bi-linear interpolation) is used in
  // the present code.
  fluid_force.assign(fluid_force.size(), {0.0, 0.0});
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (auto particle : particles) {
    for (auto &node : particle->nodes) {
      // Identify the lowest fluid lattice node in interpolation range.
      // The other fluid nodes in range have coordinates
      // (x_fluid + 1, y_fluid), (x_fluid, y_fluid + 1), and
      // (x_fluid + 1, y_fluid + 1).
      auto x_particle = node.coord[0];
      auto y_particle = node.coord[1];
      auto x_fluid = static_cast<std::size_t>(x_particle);
      auto y_fluid = static_cast<std::size_t>(y_particle);
      // Consider unrolling these 2 loops
      for (auto x = x_fluid; x < x_fluid + 2; ++x) {
        for (auto y = y_fluid; y < y_fluid + 2; ++y) {
          auto n = (y % ny) * nx + x % nx;
          // Compute interpolation weights based on distance using interpolation
          // stencil, still need to implement others
          auto weight = Dirac2(x_particle - x, y_particle - y);
          // Compute fluid force
          fluid_force[n][0] += node.force[0] * weight;
          fluid_force[n][1] += node.force[1] * weight;
        }  // y
      }  // x
    }  // node
  }  // particle
}
