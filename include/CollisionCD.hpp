#ifndef COLLISIONCD_HPP_
#define COLLISIONCD_HPP_
#include <vector>
#include "Collision.hpp"
#include "LatticeModel.hpp"

class CollisionCD: public Collision {
 public:
  /**
   * Constructor:
   * \param lat std::vector<std::vector<double>>&
   *
   */
  CollisionCD(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &position
    , const std::vector<double> &strength
    , double diffusion_coefficient
    , double initial_density_g
    , const std::vector<double> &initial_velocity);

  ~CollisionCD() = default;

  /**
   * Initializes the source lattice
   * \param position source position information
   * \param strength source magnitude at the position
   */
  void InitSource(
      const std::vector<std::vector<std::size_t>> &position
    , const std::vector<double> &strength);

  /**
   * Applies force/source term according to "A new scheme for source term in
   * LBGK model for convection-diffusion equation"
   * \param lattice 2D vector containing distribution functions
   */
  void ApplyForce(std::vector<std::vector<double>> &lattice);

  /**
   * Source term for CD equation stored row-wise
   */
  std::vector<double> source;
};
#endif  // COLLISIONCD_HPP_
