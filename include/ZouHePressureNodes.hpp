#ifndef ZOU_HE_PRESSURE_NODES_HPP_
#define ZOU_HE_PRESSURE_NODES_HPP_
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "ValueNode.hpp"

class ZouHePressureNodes: public BoundaryNodes {
 public:
  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  ZouHePressureNodes(LatticeModel &lm
    , CollisionModel &cm);

  /**
   * Destructor
   */
  ~ZouHePressureNodes() = default;

  /**
   * Assume node velocity along the side of the boundary wall is zero, i.e.,
   * u_x = 0 on top and bottom, u_y = 0 on left and right
   * \param
   * \param
   * \return
   *
   */
  void AddNode(std::size_t x
    , std::size_t y
    , double rho_node);

  void UpdateNodes(std::vector<std::vector<double>> &df
    , bool is_modify_stream);

  void UpdateSide(std::vector<std::vector<double>> &df
    , ValueNode &node);

  void UpdateCorner(std::vector<std::vector<double>> &df
    , ValueNode &node);

  /** \brief
   *
   */
  std::vector<ValueNode> nodes;

 protected:
  CollisionModel &cm_;
};
#endif  // ZOU_HE_PRESSURE_NODES_HPP_
