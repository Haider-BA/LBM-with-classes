#include "CollisionNSFCD.hpp"

CollisionNSFCD::CollisionNSFCD(LatticeModel &lm
  , const std::vector<std::vector<std::size_t>> &source_position_f
  , const std::vector<std::vector<double>> &source_strength_f
  , double kinematic_viscosity
  , double initial_density_f
  , const std::vector<std::vector<std::size_t>> &source_position_g
  , const std::vector<double> &source_strength_g
  , double diffusion_coefficient
  , double initial_density_g)
  : CollisionNSF(lm, source_position_f, source_strength_f, kinematic_viscosity,
        initial_density_f),
    CollisionCD(lm, source_position_g, source_strength_g, diffusion_coefficient,
        initial_density_g)
{}
