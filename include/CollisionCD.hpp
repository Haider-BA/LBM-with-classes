#ifndef COLLISION_CD_HPP_
#define COLLISION_CD_HPP_
#include <vector>
#include "CollisionModel.hpp"

class CollisionCD: public CollisionModel {
 public:
  /**
   * Constructor: Creates collision model for CD equation
   * \param lm lattice model used for simulation
   * \param source_position source positions
   * \param source_strength source strengths
   * \param diffusion_coefficient
   * \param initial_density_g initial density of CD lattice
   */
  CollisionCD(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<double> &source_strength
    , double diffusion_coefficient
    , double initial_density_g
    , bool is_instant);

  /**
   * Destructor
   */
  ~CollisionCD() = default;

  /**
   * Initializes the source lattice
   * \param position source position information
   * \param strength source magnitude at the position
   */
  void InitSource(
      const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<double> &source_strength);

  /**
   * Computes the macroscopic properties based on the convection-diffusion
   * collision model, just lattice density in this case. Based on "A new scheme
   * for source term in LBGK model for convection–diffusion equation
   * This is used to unify function calling in the LatticeBoltzmann TakeStep()
   * method
   * \param df lattice distribution functions stored row-wise in a 2D vector
   */
  void ComputeMacroscopicProperties(
      const std::vector<std::vector<double>> &df);

  /**
   * Applies force/source term according to "A new scheme for source term in
   * LBGK model for convection-diffusion equation"
   * \param lattice 2D vector containing distribution functions
   */
  void Collide(std::vector<std::vector<double>> &df);

  /**
   * Sets source term to 0
   */
  void KillSource();

  /**
   * Source term for CD equation stored row-wise
   */
  std::vector<double> source;

 protected:
  /**
   * Boolean toggle to indicate if the source term in this collision model is an
   * instantaneous source, for use with diffusion analytical solution
   */
  bool is_instant_;
};
#endif  // COLLISION_CD_HPP_
