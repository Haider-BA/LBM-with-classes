#include "CollisionCD.hpp"

CollisionCD::CollisionCD(double initial_density_g)
  : rho0_ {initial_density_g}
{
  is_implemented = true;
}
