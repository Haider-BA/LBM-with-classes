#ifndef COLLISIONNS_HPP_
#define COLLISIONNS_HPP_
#include <vector>
#include "Collision.hpp"
#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"

class CollisionNS: public Collision {
 public:
  /**
   * Constructor: Creates collision model for NS equation
   * \param lm lattice model used for simulation
   * \param position source positions
   * \param strenght source strengths
   * \param kinematic viscosity
   * \param initial_density_g initial density of NS lattice
   */
  CollisionNS(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &position
    , const std::vector<std::vector<double>> &strength
    , double kinematic_viscosity
    , double initial_density_f);

  /**
   * Destructor
   */
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
  void Collide(std::vector<std::vector<double>> &lattice);

  /**
   * Sets source term to 0
   */
  void KillSource();

  /**
   * Source term for NS equation stored row-wise
   */
  std::vector<std::vector<double>> source;
};
#endif  // COLLISIONNS_HPP_
