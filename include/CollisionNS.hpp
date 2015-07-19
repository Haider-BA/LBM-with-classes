#ifndef COLLISION_NS_HPP_
#define COLLISION_NS_HPP_
#include <vector>
#include "CollisionModel.hpp"
#include "LatticeModel.hpp"

class CollisionNS: public CollisionModel {
 public:
  /**
   * Constructor: Creates collision model for NS equation with the same density
   * at each node
   * \param lm lattice model used for simulation
   * \param kinematic viscosity
   * \param initial_density_f initial density of NS lattice
   */
  CollisionNS(LatticeModel &lm
    , double kinematic_viscosity
    , double initial_density_f);

  /**
   * Destructor
   */
  virtual ~CollisionNS() = default;

  /**
   * Collides according to Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  virtual void Collide(std::vector<std::vector<double>> &lattice);
};

#endif  // COLLISION_NS_HPP_
