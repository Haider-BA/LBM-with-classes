#ifndef COLLISIONNS_HPP_
#define COLLISIONNS_HPP_
#include <vector>
#include "Collision.hpp"
#include "LatticeModel.hpp"

class CollisionNS: public Collision {
 public:
  /** \brief
   *
   * \param lat std::vector<std::vector<double>>&
   *
   */
  CollisionNS(LatticeModel &lm
    , const std::vector<std::vector<std::size_t>> &position
    , const std::vector<std::vector<double>> &strength
    , double kinematic_viscosity
    , double initial_density_f
    , const std::vector<double> &initial_velocity);

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  void InitSource(
      const std::vector<std::vector<std::size_t>> &position
    , const std::vector<std::vector<double>> &strength);

  /** \brief
   *
   * \return void
   *
   */
  void ApplyForce();
  // only public for ease of debugging with Print() function
  std::vector<std::vector<double>> source_;
 private:
//  std::vector<std::vector<double>> source_;
};
#endif  // COLLISIONNS_HPP_
