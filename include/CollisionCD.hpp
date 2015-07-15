#ifndef COLLISIONCD_HPP_
#define COLLISIONCD_HPP_
#include "Collision.hpp"

class CollisionCD: public Collision{
 public:
  CollisionCD() = default;
  CollisionCD(double initial_density_g);
  double rho0_;
 private:
  // double rho0_;
};
#endif  // COLLISIONCD_HPP_
