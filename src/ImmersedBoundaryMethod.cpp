#include "ImmersedBoundaryMethod.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <vector>
#include "Algorithm.hpp"
#include "LatticeModel.hpp"

ImmersedBoundaryMethod::ImmersedBoundaryMethod(int stencil
  , std::vector<std::vector<double>> &lattice_force
  , LatticeModel &lm)
  : particles {},
    fluid_force (lattice_force),
    stencil_ {stencil},
    lm_ (lm)
{}

void ImmersedBoundaryMethod::AddParticle(Particle *particle)
{
  particles.push_back(particle);
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
      const auto x_particle = node.coord[0];
      const auto y_particle = node.coord[1];
      const auto x_fluid = static_cast<int>(floor(x_particle + 2 *
          particle->char_length * std::numeric_limits<double>::epsilon()));
      const auto y_fluid = static_cast<int>(floor(y_particle + 2 *
          particle->char_length * std::numeric_limits<double>::epsilon()));
      std::cout << x_particle << " " << x_fluid << std::endl;
      // Consider unrolling these 2 loops
      for (auto y = y_fluid; y < y_fluid + 2; ++y) {
        for (auto x = x_fluid; x < x_fluid + 2; ++x) {
          const auto n = ((y + ny) % ny) * nx + (x + nx) % nx;
          // Compute interpolation weights based on distance using interpolation
          // stencil, still need to implement others
          const auto weight = Dirac(stencil_, x_particle - x, y_particle - y);
          // Compute fluid force
          fluid_force[n][0] += node.force[0] * scaling * weight;
          fluid_force[n][1] += node.force[1] * scaling * weight;
        }  // x
      }  // y
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
      const auto x_particle = node.coord[0];
      const auto y_particle = node.coord[1];
      // fix for truncation error when double is really close to integer value
      // http://blog.frama-c.com/index.php?post/2013/05/02/nearbyintf1
      // uses "machine epsilon" through std::numeric_limits<double>::epsilon()
      // http://stackoverflow.com/questions/17333/most-effective-way-for-float-
      // and-double-comparison
      const auto x_fluid = static_cast<int>(floor(x_particle + 2 *
          particle->char_length * std::numeric_limits<double>::epsilon()));
      const auto y_fluid = static_cast<int>(floor(y_particle + 2 *
          particle->char_length * std::numeric_limits<double>::epsilon()));
      for (auto y = y_fluid; y < y_fluid + 2; ++y) {
        for (auto x = x_fluid; x < x_fluid + 2; ++x) {
          const auto n = ((y + ny) % ny) * nx + (x + nx) % nx;
          const auto weight = Dirac(stencil_, x_particle - x, y_particle - y);
          // node velocity maybe calculated in dimensionless units (check)
          node.u[0] += lm_.u[n][0] / scaling * weight;
          node.u[1] += lm_.u[n][1] / scaling * weight;
        }  // x
      }  // y
    }  // node
  }  // particle
}

void ImmersedBoundaryMethod::UpdateParticlePosition()
{
  const auto dt = lm_.GetTimeStep();
  const auto nx = lm_.GetNumberOfColumns();
  for (auto particle : particles) {
    const auto nn = particle->GetNumberOfNodes();
    particle->center.coord = {0.0, 0.0};
    for (auto &node : particle->nodes) {
      node.coord[0] += node.u[0] * dt;
      node.coord[1] += node.u[1] * dt;
      particle->center.coord[0] += node.coord[0] / nn;
      particle->center.coord[1] += node.coord[1] / nn;
    }  // node
    if (particle->is_mobile) particle->UpdateReferencePosition();
    if (particle->center.coord[0] < 0) {
      particle->center.coord[0] += nx;
      for (auto &node : particle->nodes) node.coord[0] += nx;
      if (particle->is_mobile) {
        particle->center.coord_ref[0] += nx;
        for (auto &node : particle->nodes) node.coord_ref[0] += nx;
      }
    }
    else if (particle->center.coord[0] >= nx) {
      particle->center.coord[0] -= nx;
      for (auto &node : particle->nodes) node.coord[0] -= nx;
      if (particle->is_mobile) {
        particle->center.coord_ref[0] -= nx;
        for (auto &node : particle->nodes) node.coord_ref[0] -= nx;
      }
    }
  }  // particle
}
