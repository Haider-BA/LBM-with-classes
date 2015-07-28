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

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  void AddNode(std::size_t x
    , std::size_t y
    , double rho_node
    , double u_x
    , double v_y);

  /** \brief
   *
   */
  std::vector<ValueNode> nodes;

 protected:
  CollisionModel &cm_;
};
#endif  // ZOU_HE_PRESSURE_NODES_HPP_
