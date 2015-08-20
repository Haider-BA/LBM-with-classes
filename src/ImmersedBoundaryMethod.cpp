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
  lm_.u.assign(lm_.u.size(), {0.0, 0.0});
}

void ImmersedBoundaryMethod::SpreadForce()
{
  // The two-point interpolation stencil (bi-linear interpolation) is used in
  // the present code.
  fluid_force.assign(fluid_force.size(), {0.0, 0.0});
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (auto particle : particles) {
    auto nn = particle->GetNumberOfNodes();
    for (auto i = 0u; i < nn; ++i) {
      // Identify the lowest fluid lattice node in interpolation range.
      // The other fluid nodes in range have coordinates
      // (x_fluid + 1, y_fluid), (x_fluid, y_fluid + 1), and
      // (x_fluid + 1, y_fluid + 1).
      auto x_particle = particle->nodes[i].coord[0];
      auto y_particle = particle->nodes[i].coord[1];
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
          fluid_force[n][0] += particle->nodes[i].force[0] * weight;
          fluid_force[n][1] += particle->nodes[i].force[1] * weight;
        }  // y
      }  // x
    }  // i
  }  // particle
}
