#ifndef VALUE_NODE_HPP_
#define VALUE_NODE_HPP_
#include "Node.hpp"

class ValueNode: public Node {
 public:
  /**
   * Constructor: Creates a node which contains information on its position,
   * 2 double values, 1 boolean value and 1 integer value. To be used by
   * Zou/He velocity node
   * \param x_position x-coordinate of the node
   * \param y_position y-coordinate of the node
   * \param nx number of columns of the lattice model
   * \param d_1 first double value to be stored in vector of double values, used
   *        as x-velocity of Zou/He velocity node
   * \param d_2 second double value to be stored in vector of double values,
   *        used as y-velocity of Zou/He velocity node
   * \param b_1 first boolean value, used to indicate if node is a corner node
   *        in Zou/He velocity node
   * \param i_1 first integer value, used to indicate which side the node
   *        belongs to for non-corner nodes, and which corner the node belongs
   *        to for corner nodes in Zou/He velocity nodes
   */
  ValueNode(std::size_t x_position
    , std::size_t y_position
    , std::size_t nx
    , double d_1
    , double d_2
    , bool b_1
    , int i_1);

  /**
   * Constructor: Creates a node which contains information on its position,
   * 1 double value, 1 boolean value and 1 integer value. To be used by
   * Zou/He pressure node
   * \param x_position x-coordinate of the node
   * \param y_position y-coordinate of the node
   * \param nx number of columns of the lattice model
   * \param d_1 first double value, used as node pressure in Zou/He pressure
   *        nodes
   * \param b_1 first boolean value, used to indicate if node is a corner node
   *        Zou/He pressure nodes
   * \param i_1 first integer value, used to indicate which side the node
   *        belongs to for non-corner nodes, and which corner the node belongs
   *        to for corner nodes in Zou/He pressure nodes
   */
  ValueNode(std::size_t x_position
    , std::size_t y_position
    , std::size_t nx
    , double d_1
    , bool b_1
    , int i_1);

  /**
   * Destructor
   */
  ~ValueNode() = default;

  /**
   * Double value, used as node pressure in Zou/He pressure nodes
   */
  double d1;

  /**
   * Vector containing double values, used as node velocity in Zou/He velocity
   * nodes
   */
  std::vector<double> v1;

  /**
   * Boolean value, used to indicate if node is a corner node Zou/He velocity
   * and pressure nodes
   */
  bool b1;

  /**
   * Integer value, used to indicate which side the node belongs to for
   * non-corner nodes, and which corner the node belongs to for corner nodes in
   * Zou/He velocity and pressure nodes
   */
  int i1;
};

#endif  // VALUE_NODE_HPP_
