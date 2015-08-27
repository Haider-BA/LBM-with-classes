#ifndef PARTICLE_NODE_HPP_
#define PARTICLE_NODE_HPP_
#include <vector>

class ParticleNode {
 public:
  /**
   * Particle node without velocity and force components, for use as particle
   * centers
   * \param x x-coordinate
   * \param y y-coordinate
   * \param x_ref reference x_coordinate
   * \param y_ref reference y_coordinate
   */
  ParticleNode(double x
    , double y
    , double x_ref
    , double y_ref);

  /**
   * Particle node with velocity and force components, for use as regular
   * particle nodes
   * \param x x-coordinate
   * \param y y-coordinate
   * \param x_ref reference x_coordinate
   * \param y_ref reference y_coordinate
   * \param u_x velocity in x-direction
   * \param u_y velocity in y-direction
   * \param force_x force in x-direction
   * \param force_y force in y-direction
   */
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
