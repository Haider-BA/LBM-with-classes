#include "ImmersedBoundaryMethod.hpp"

ImmersedBoundaryMethod::ImmersedBoundaryMethod(int interpolation_stencil
  , std::vector<std::vector<double>> &lattice_force
  , std::vector<std::vector<double>> &lattice_velocity)
  : particles{},
    fluid_force {lattice_force},
    fluid_velocity {lattice_velocity},
    interpolation_stencil_ {interpolation_stencil}
{}

void ImmersedBoundaryMethod::InterpolateFluidVelocity()
{
  fluid_velocity.assign(fluid_velocity.size(), {0.0, 0.0});
}

void ImmersedBoundaryMethod::SpreadForce()
{
  fluid_force.assign(fluid_force.size(), {0.0, 0.0});
}
