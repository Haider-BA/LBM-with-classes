#ifndef ZOU_HE_PRESSURE_NODES_HPP_
#define ZOU_HE_PRESSURE_NODES_HPP_
#include <vector>
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "ValueNode.hpp"

class ZouHePressureNodes: public BoundaryNodes {
 public:
  /**
   * Constructor: Creates Zou/He pressure boundary nodes
   * \param lm lattice model which contains information on the number of rows,
   *        columns, dimensions, discrete directions and lattice velocity
   * \param cm collision model which contains information on lattice density
   */
  ZouHePressureNodes(LatticeModel &lm
    , CollisionModel &cm);

  /**
   * Destructor
   */
  ~ZouHePressureNodes() = default;

  /**
   * Adds a Zou/He pressure node to the nodes vector. Assume node velocity along
   * the side of the boundary wall is zero, i.e., u_x = 0 on top and bottom,
   * u_y = 0 on left and right
   * \param x x-coordinate of the node
   * \param y y-coordinate of the node
   * \param rho_node pressure/density of the node
   */
  void AddNode(std::size_t x
    , std::size_t y
    , double rho_node);

  /**
   * Updates the boundary nodes based on "On pressure and velocity boundary
   * conditions for the lattice Boltzmann"
   * \param df lattice distribution functions stored row-wise in a 2D vector
   * \param is_modify_stream boolean toggle for half-way bounceback nodes to
   *        perform functions during stream, set to FALSE for Zou/He pressure
   *        nodes
   */
  void UpdateNodes(std::vector<std::vector<double>> &df
    , bool is_modify_stream);

  /**
   * Updates the non-corner nodes
   * \param df lattice distribution function stored row-wise in a 2D vector
   * \param node Zou/He pressure node which contains information on the position
   *        of the boundary node, pressure/density and velocities of the node
   */
  void UpdateSide(std::vector<std::vector<double>> &df
    , ValueNode &node);

  /**
   * Updates the corner nodes, first-order expolation for node density
   * \param df lattice distribution function stored row-wise in a 2D vector
   * \param node Zou/He pressure node which contains information on the position
   *        of the boundary node, pressure/density and velocities of the node
   */
  void UpdateCorner(std::vector<std::vector<double>> &df
    , ValueNode &node);

  /**
   * Boundary nodes stored in a 1D vector
   */
  std::vector<ValueNode> nodes;

 protected:
  /**
   * Collision model which contains information on lattice density
   */
  CollisionModel &cm_;

  /**
   * Additional constants beta1, beta2 and beta3 since the dx = dt = 1 condition
   * is not always maintained
   */
  double beta1_;
  double beta2_;
  double beta3_;
};
#endif  // ZOU_HE_PRESSURE_NODES_HPP_
