#ifndef ZOU_HE_NODES_HPP_
#define ZOU_HE_NODES_HPP_
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "ValueNode.hpp"

class ZouHeNodes: public BoundaryNodes {
 public:
  /** \brief
   *
   * \param
   * \param
   */
  ZouHeNodes(LatticeModel &lm
    , CollisionModel &cm);

  /**
   * Destructor
   */
  ~ZouHeNodes() = default;

  /** \brief
   *
   * \param
   * \param
   * \param
   * \param
   *
   */
  void AddNode(std::size_t x
    , std::size_t y
    , double u_x
    , double u_y);

  /** \brief
   *
   * \param
   * \param
   */
  void UpdateNodes(std::vector<std::vector<double>> &df
    , bool is_modify_stream);

  /** \brief
   *
   * \param
   * \param
   */
  void UpdateSide(std::vector<std::vector<double>> &df
    , ValueNode &node);

  /** \brief
   *
   * \param
   * \param
   */
  void UpdateCorner(std::vector<std::vector<double>> &df
    , ValueNode &node);

  /** \brief
   *
   */
  std::vector<ValueNode> nodes;

 protected:
  /** \brief
   *
   */
  CollisionModel &cm_;
//  std::vector<bool> is_corner_;
//  std::vector<bool> knowns_;

};

#endif  // ZOU_HE_NODES_HPP_
