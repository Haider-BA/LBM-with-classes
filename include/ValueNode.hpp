#ifndef VALUE_NODE_HPP_
#define VALUE_NODE_HPP_
#include "Node.hpp"

class ValueNode: public Node {
 public:
  ValueNode(std::size_t x_position
    , std::size_t y_position
    , std::size_t nx
    , double d_1
    , double d_2
    , bool b_1
    , int i_1);

  ValueNode(std::size_t x_position
    , std::size_t y_position
    , std::size_t nx
    , double d_1
    , bool b_1
    , int i_1);

  ~ValueNode() = default;

  double d1;
  std::vector<double> v1;
  bool b1;
  int i1;
};

#endif  // VALUE_NODE_HPP_
