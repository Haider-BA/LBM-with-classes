#ifndef COLLISIONCDNS_HPP_
#define COLIISIONCDNS_HPP_
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"

class CollisionCDNS: public CollisionCD, public CollisionNS {
 public:
  CollisionCDNS(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<double> &source_strength
    , double diffusion_coefficient
    , double initial_density_g
    , double kinematic_viscosity
    , double initial_density_f);

  ~CollisionCDNS() = default;
};

#endif  // COLLISIONCDNS_HPP_
