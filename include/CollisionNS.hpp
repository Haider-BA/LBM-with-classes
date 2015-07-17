#ifndef COLLISIONNS_HPP_
#define COLLISIONNS_HPP_
#include <vector>
#include "CollisionModel.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"

class CollisionNS: public CollisionModel {
 public:
  /**
   * Constructor: Creates collision model for NS equation with the same density
   * at each node
   * \param lm lattice model used for simulation
   * \param position source positions
   * \param strenght source strengths
   * \param kinematic viscosity
   * \param initial_density_f initial density of NS lattice
   */
  CollisionNS(LatticeModel &lm
    , double kinematic_viscosity
    , double initial_density_f);

  /**
   * Constructor: Creates collision model for NS equation with variable density
   * at each node
   * \param lm lattice model used for simulation
   * \param position source positions
   * \param strenght source strengths
   * \param kinematic viscosity
   * \param initial_density_f initial density of NS lattice
   */
  CollisionNS(LatticeModel &lm
    , double kinematic_viscosity
    , const std::vector<double> &initial_density_f);

  /**
   * Destructor
   */
  ~CollisionNS() = default;

  /**
   * Collides according to Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  virtual void CollideNS(std::vector<std::vector<double>> &lattice);

  /**
   * Collision for Lid-driven flow, to be refactored
   */
  void CollideLid(std::vector<std::vector<double>> &lattice);

 protected:
  double tau_f_;
};
#endif  // COLLISIONNS_HPP_
