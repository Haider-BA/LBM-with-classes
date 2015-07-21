#ifndef COLLISION_NSF_HPP_
#define COLLISION_NSF_HPP_
#include <vector>
#include "CollisionNS.hpp"

class CollisionNSF: public CollisionNS {
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
  CollisionNSF(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<std::vector<double>> &source_strength
    , double kinematic_viscosity
    , double initial_density_f);

  /**
   * Destructor
   */
  ~CollisionNSF() = default;

  /**
   * Initializes the source lattice
   * \param position source position information
   * \param strength source magnitude at the position in the various dimensions
   */
  void InitSource(
      const std::vector<std::vector<std::size_t>> &source_position
    , const std::vector<std::vector<double>> &source_strength);

  /**
   * Calculated velocity for NS equation based on formula in
   * Guo2002
   * \param df 2D vector containing distribution functions of the NS
   *        equation
   * \return 2D vector containing velocity at each node of lattice
   */
  std::vector<std::vector<double>> ComputeU(
      const std::vector<std::vector<double>> &df);

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  void ComputeMacroscopicProperties(
      const std::vector<std::vector<double>> &df);

  /**
   * Collides and applies force according to Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  void Collide(std::vector<std::vector<double>> &lattice);

  /**
   * Source term for NS equation stored row-wise
   */
  std::vector<std::vector<double>> source;
};
#endif // COLLISION_NSF_HPP_
