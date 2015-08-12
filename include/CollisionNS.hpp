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
   * Constructor: Creates collision model for NS equation with the same density
   * at each node
   * \param lm lattice model used for simulation
   * \param kinematic viscosity
   * \param initial_density_f initial density of NS lattice
   */
  CollisionNS(LatticeModel &lm
    , double kinematic_viscosity
    , const std::vector<double> &initial_density_f);

  /**
   * Virtual destructor since we may be deriving from this class
   */
  virtual ~CollisionNS() = default;

  /**
   * Calculated velocity for NS equation without body force based on formula in
   * Guo2002
   * \param df 2D vector containing distribution functions of the NS
   *        equation
   * \return 2D vector containing velocity at each node of lattice
   */
  virtual std::vector<std::vector<double>> ComputeU(
      const std::vector<std::vector<double>> &df);

  /**
   * Computes the macroscopic properties based on the collision model used, both
   * velocity and density in this case. Based on "Discrete lattice effects on
   * the forcing term in the lattice Boltzmann method"
   * \param df lattice distribution functions stored row-wise in a 2D vector
   */
  void ComputeMacroscopicProperties(
      const std::vector<std::vector<double>> &df);

  /**
   * Collides according to Guo2002
   * \param lattice 2D vector containing distribution functions
   */
  virtual void Collide(std::vector<std::vector<double>> &lattice);
};

#endif  // COLLISION_NS_HPP_
