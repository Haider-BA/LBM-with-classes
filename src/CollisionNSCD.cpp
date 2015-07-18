#include "CollisionNSCD.hpp"

CollisionNSCD::CollisionNSCD(LatticeModel &lm
  , double kinematic_viscosity
  , double initial_density_f
  , const std::vector<std::vector<std::size_t>> &source_position
  , const std::vector<double> &source_strength
  , double diffusion_coefficient
  , double initial_density_g)
  : CollisionNS(lm, kinematic_viscosity, initial_density_f),
    CollisionCD(lm, source_position, source_strength, diffusion_coefficient,
        initial_density_g)
{}
