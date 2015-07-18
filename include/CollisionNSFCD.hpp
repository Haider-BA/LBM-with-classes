#ifndef COLLISIONNSFCD_HPP_
#define COLLISIONNSFCD_HPP_
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNSF.hpp"

class CollisionNSFCD: public CollisionNSF, public CollisionCD {
 public:
  CollisionNSFCD(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &source_position_f
    , const std::vector<std::vector<double>> &source_strength_f
    , double kinematic_viscosity
    , double initial_density_f
    , const std::vector<std::vector<std::size_t>> &source_position_g
    , const std::vector<double> &source_strength_g
    , double diffusion_coefficient
    , double initial_density_g);

  ~CollisionNSFCD() = default;
};

#endif  // COLLISIONNSFCD_HPP_
