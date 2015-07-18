#ifndef COLLISIONNSCD_HPP_
#define COLIISIONNSCD_HPP_
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"

class CollisionNSCD: public CollisionNS, public CollisionCD {
 public:
  CollisionNSCD(LatticeModel &lm
    , double kinematic_viscosity
    , double initial_density_f
    , const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<double> &source_strength
    , double diffusion_coefficient
    , double initial_density_g);

  ~CollisionNSCD() = default;
};

#endif  // COLLISIONNSCD_HPP_
