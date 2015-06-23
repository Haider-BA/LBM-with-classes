#ifndef COLLISIONNS_HPP_
#define COLLISIONNS_HPP_
#include <vector>
#include "Collision.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"

class CollisionNS: public Collision {
 public:
  /**
   * Constructor:
   * \param lat std::vector<std::vector<double>>&
   *
   */
  CollisionNS(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &position
    , const std::vector<std::vector<double>> &strength
    , double kinematic_viscosity
    , double initial_density_f
    , const std::vector<double> &initial_velocity);

  ~CollisionNS() = default;

  /**
   * Initializes the source lattice
   * \param position source position information
   * \param strength source magnitude at the position in the various dimensions
   */
  void InitSource(
      const std::vector<std::vector<std::size_t>> &position
    , const std::vector<std::vector<double>> &strength);

  /**
   * Applies force according to Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  void ApplyForce(std::vector<std::vector<double>> &lattice);

  /**
   * Source term for NS equation stored row-wise
   */
  std::vector<std::vector<double>> source;
};
#endif  // COLLISIONNS_HPP_
