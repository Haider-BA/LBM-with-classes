#ifndef COLLISIONCD_HPP_
#define COLLISIONCD_HPP_
#include <vector>
#include "CollisionModel.hpp"
#include "LatticeModel.hpp"

class CollisionCD: public CollisionModel {
 public:
  /**
   * Constructor: Creates collision model for CD equation
   * \param lm lattice model used for simulation
   * \param position source positions
   * \param strenght source strengths
   * \param diffusion_coefficient
   * \param initial_density_g initial density of CD lattice
   */
  CollisionCD(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<double> &source_strength
    , double diffusion_coefficient
    , double initial_density_g);

  /**
   * Destructor
   */
  virtual ~CollisionCD() = default;

  /**
   * Initializes the source lattice
   * \param position source position information
   * \param strength source magnitude at the position
   */
  void InitSource(
      const std::vector<std::vector<std::size_t>> &position
    , const std::vector<double> &strength);

  /**
   * Collides and applies force/source term according to "A new scheme for
   * source term in LBGK model for convection-diffusion equation"
   * \param lattice 2D vector containing distribution functions
   */
  void CollideCD(std::vector<std::vector<double>> &lattice);

  /**
   * Sets source term to 0
   */
  void KillSource();

  /**
   * Density stored row-wise in a 1D vector.
   */
  std::vector<double> rho_g;

  /**
   * Source term for CD equation stored row-wise
   */
  std::vector<double> source;

 protected:
  double tau_g_;
};
#endif  // COLLISIONCD_HPP_
