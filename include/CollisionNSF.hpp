#ifndef COLLISIONNSF_HPP_
#define COLLISIONNSF_HPP_

#include "CollisionNS.hpp"

class CollisionNSF: public CollisionNS {
 public:
  CollisionNSF(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<std::vector<double>> &source_strength
    , double kinematic_viscosity
    , double initial_density_f);

  /**
   * Initializes the source lattice
   * \param position source position information
   * \param strength source magnitude at the position in the various dimensions
   */
  void InitSource(
      const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<std::vector<double>> &source_strength);

  /**
   * Collides and applies force according to Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  void CollideNS(std::vector<std::vector<double>> &lattice);

  /**
   * Sets source term to 0
   */
  void KillSource();

  /**
   * Source term for NS equation stored row-wise
   */
  std::vector<std::vector<double>> source;
};
#endif // COLLISIONNSF_HPP_
