#include "ImmersedBoundaryMethod.hpp"
#include <iostream>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

ImmersedBoundaryMethod::ImmersedBoundaryMethod(int interpolation_stencil
  , std::vector<std::vector<double>> &lattice_force
  , LatticeModel &lm)
  : particles {},
    fluid_force (lattice_force),
    interpolation_stencil_ {interpolation_stencil},
    lm_ (lm)
{}

void ImmersedBoundaryMethod::AddParticle(Particle *particle)
{
  particles.push_back(particle);
}

double ImmersedBoundaryMethod::Phi2(double x
  , double h)
{
  double phi = 0;
  auto x_abs = fabs(x);
  if (x_abs <= h) phi = (1 - x_abs / h);
  return phi;
}

double ImmersedBoundaryMethod::Phi3(double x
  , double h)
{
  double phi = 0;
  auto x_abs = fabs(x);
  if (x_abs <= 0.5) {
    phi = (1 + sqrt(1 - 3 * x * x)) / 3;
  }
  else if (x_abs <= 1.5) {
    phi = (5 - 3 * x_abs - sqrt(-2 + 6 * x_abs - 3 * x * x)) / 6;
  }
  return phi;
}

double ImmersedBoundaryMethod::Phi4(double x
  , double h)
{
  double phi = 0;
  auto x_abs = fabs(x);
  if (x_abs <= 1) {
    phi = (3 - 2 * x_abs + sqrt(1 + 4 * x_abs - 4 * x * x)) / 8;
  }
  else if (x_abs <= 2) {
    phi = (5 - 2 * x_abs - sqrt(-7 + 12 * x_abs - 4 * x * x)) / 8;
  }
  return phi;
}

double ImmersedBoundaryMethod::Dirac2(double x
  , double y
  , double h)
{
  return ImmersedBoundaryMethod::Phi2(x, h) *
      ImmersedBoundaryMethod::Phi2(y, h);
}

double ImmersedBoundaryMethod::Dirac3(double x
  , double y
  , double h)
{
  return ImmersedBoundaryMethod::Phi3(x, h) *
      ImmersedBoundaryMethod::Phi3(y, h);
}

double ImmersedBoundaryMethod::Dirac4(double x
  , double y
  , double h)
{
  return ImmersedBoundaryMethod::Phi4(x, h) *
      ImmersedBoundaryMethod::Phi4(y, h);
}

void ImmersedBoundaryMethod::SpreadForce()
{
  // The two-point interpolation stencil (bi-linear interpolation) is used in
  // the present code.
  fluid_force.assign(fluid_force.size(), {0.0, 0.0});
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  const auto dx = lm_.GetSpaceStep();
  const auto dt = lm_.GetTimeStep();
  const auto scaling = dx / dt / dt;
  for (auto particle : particles) {
    for (auto &node : particle->nodes) {
      // Identify the lowest fluid lattice node in interpolation range.
      // The other fluid nodes in range have coordinates
      // (x_fluid + 1, y_fluid), (x_fluid, y_fluid + 1), and
      // (x_fluid + 1, y_fluid + 1).
      const auto x_particle = node.coord[0] / dx;
      const auto y_particle = node.coord[1] / dx;
      const auto x_fluid = static_cast<std::size_t>(x_particle);
      const auto y_fluid = static_cast<std::size_t>(y_particle);
      // Consider unrolling these 2 loops
      for (auto x = x_fluid; x < x_fluid + 2; ++x) {
        for (auto y = y_fluid; y < y_fluid + 2; ++y) {
          const auto n = (y % ny) * nx + x % nx;
          // Compute interpolation weights based on distance using interpolation
          // stencil, still need to implement others
          const auto weight = ImmersedBoundaryMethod::Dirac2((x_particle - x) *
              dx, (y_particle - y) * dx, dx);
          // Compute fluid force
          fluid_force[n][0] += node.force[0] * scaling * weight;
          fluid_force[n][1] += node.force[1] * scaling * weight;
        }  // y
      }  // x
    }  // node
  }  // particle
}

void ImmersedBoundaryMethod::InterpolateFluidVelocity()
{
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  const auto dx = lm_.GetSpaceStep();
  const auto dt = lm_.GetTimeStep();
  const auto scaling = dx / dt;
  // pointer always editable, do not need &
  for (auto particle : particles) {
    for (auto &node : particle->nodes) {
      node.u = {0.0, 0.0};
      const auto x_particle = node.coord[0] / dx;
      const auto y_particle = node.coord[1] / dx;
      const auto x_fluid = static_cast<std::size_t>(x_particle);
      const auto y_fluid = static_cast<std::size_t>(y_particle);
      for (auto x = x_fluid; x < x_fluid + 2; ++x) {
        for (auto y = y_fluid; y < y_fluid + 2; ++y) {
          const auto n = (y % ny) * nx + x % nx;
          const auto weight = ImmersedBoundaryMethod::Dirac2((x_particle - x) *
              dx, (y_particle - y) * dx, dx);
          // node velocity maybe calculated in dimensionless units (check)
          node.u[0] += lm_.u[n][0] / scaling * weight;
          node.u[1] += lm_.u[n][1] / scaling * weight;
        }  // y
      }  // x
    }  // node
  }  // particle
}

void ImmersedBoundaryMethod::UpdateParticlePosition()
{
  const auto dt = lm_.GetTimeStep();
  for (auto particle : particles) {
    const auto nn = particle->GetNumberOfNodes();
    particle->center.coord = {0.0, 0.0};
    for (auto &node : particle->nodes) {
      node.coord[0] += node.u[0] * dt;
      node.coord[1] += node.u[1] * dt;
      particle->center.coord[0] += node.coord[0] / nn;
      particle->center.coord[1] += node.coord[1] / nn;
    }  // node
  }  // particle
}
