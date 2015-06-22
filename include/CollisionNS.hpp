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

  ~CollisionNS() = default;

  // maybe use friend class to access private variables easier for testing?
  std::vector<std::vector<double>> GetSource() const;

  /**
   * Initializes the source lattice
   * \param position source position information
   * \param strength source magnitude at the position in the various dimensions
   */
  void InitSource(
      const std::vector<std::vector<std::size_t>> &position
    , const std::vector<std::vector<double>> &strength);

  /** \brief
   *
   * \return void
   *
   */
  void ApplyForce(std::vector<std::vector<double>> &lattice);

 private:
  std::vector<std::vector<double>> source_;
};
#endif  // COLLISIONNS_HPP_
