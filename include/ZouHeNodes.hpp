#ifndef ZOU_HE_NODES_HPP_
#define ZOU_HE_NODES_HPP_
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "ValueNode.hpp"

class ZouHeNodes: public BoundaryNodes {
 public:
  /**
   * Constructor: Creates Zou/He velocity boundary nodes
   * \param lm lattice model which contains information on the number of rows,
   *        columns, dimensions, discrete directions and lattice velocity
   * \param cm collision model which contains information on lattice density
   */
  ZouHeNodes(LatticeModel &lm
    , CollisionModel &cm);

  /**
   * Destructor
   */
  ~ZouHeNodes() = default;

  /**
   * Adds a Zou/He velocity node to the nodes vector
   * \param x x-coordinate of the node
   * \param y y-coordinate of the node
   * \param u_x x-velocity of the node
   * \param u_y y-velocity of the node
   */
  void AddNode(std::size_t x
    , std::size_t y
    , double u_x
    , double u_y);

  /**
   * Updates the boundary nodes based on "On pressure and velocity boundary
   * conditions for the lattice Boltzmann"
   * \param df lattice distribution functions stored row-wise in a 2D vector
   * \param is_modify_stream boolean toggle for half-way bounceback nodes to
   *        perform functions during stream, set to FALSE for Zou/He velocity
   *        nodes
   */
  void UpdateNodes(std::vector<std::vector<double>> &df
    , bool is_modify_stream);

  /**
   * Updates the non-corner nodes
   * \param df lattice distribution function stored row-wise in a 2D vector
   * \param node Zou/He velocity node which contains information on the position
   *        of the boundary node and velocities of the node
   */
  void UpdateSide(std::vector<std::vector<double>> &df
    , ValueNode &node);

  /**
   * Updates the corner nodes, first-order expolation for node density
   * \param df lattice distribution function stored row-wise in a 2D vector
   * \param node Zou/He velocity node which contains information on the position
   *        of the boundary node and velocities of the node
   */
  void UpdateCorner(std::vector<std::vector<double>> &df
    , ValueNode &node);

  /**
   * Toggles behaviour of Zou/He nodes when used as outlet, boundary node
   * velocity will be extrapolated (1st order) from the neighbouring nodes
   */
  void ToggleNormalFlow();

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
   * Boolean toggle for open boundary condition (outlet)
   */
  bool is_normal_flow_;

  /**
   * Additional constants beta1, beta2 and beta3 since the dx = dt = 1 condition
   * is not always maintained
   */
  double beta1_;
  double beta2_;
  double beta3_;
//  std::vector<bool> is_corner_;
//  std::vector<bool> knowns_;

};

#endif  // ZOU_HE_NODES_HPP_
