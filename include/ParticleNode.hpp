#ifndef PARTICLE_NODE_HPP_
#define PARTICLE_NODE_HPP_
#include <vector>

class ParticleNode {
 public:
  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  ParticleNode(double x
    , double y
    , double x_ref
    , double y_ref);

  ParticleNode(double x
    , double y
    , double x_ref
    , double y_ref
    , double u_x
    , double u_y
    , double force_x
    , double force_y);

  ~ParticleNode() = default;

  std::vector<double> coord;
  std::vector<double> coord_ref;
  std::vector<double> u;
  std::vector<double> force;
};
#endif  // PARTICLE_NODE_HPP_
