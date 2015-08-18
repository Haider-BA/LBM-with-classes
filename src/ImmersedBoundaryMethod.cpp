#include "ImmersedBoundaryMethod.hpp"
#include <iostream>

ImmersedBoundaryMethod::ImmersedBoundaryMethod(int interpolation_stencil
  , std::vector<std::vector<double>> &lattice_force
  , std::vector<std::vector<double>> &lattice_velocity)
  : particles{},
    fluid_force {lattice_force},
    fluid_velocity {lattice_velocity},
    interpolation_stencil_ {interpolation_stencil}
{}

void ImmersedBoundaryMethod::AddParticle(Particle* particle)
{
  particles.push_back(particle);
}

void ImmersedBoundaryMethod::InterpolateFluidVelocity()
{
  fluid_velocity.assign(fluid_velocity.size(), {0.0, 0.0});
}

void ImmersedBoundaryMethod::SpreadForce()
{
  fluid_force.assign(fluid_force.size(), {0.0, 0.0});
  for (auto particle : particles) {
    std::cout << "looping particles" << std::endl;
    auto nn = particle->GetNumberOfNodes();
    for (auto n = 0u; n < nn; ++n) {
      std::cout << "node " << n << " "
                << particle->nodes[n].coord[0] << " "
                << particle->nodes[n].coord[1] << std::endl;
    }
  }
}
