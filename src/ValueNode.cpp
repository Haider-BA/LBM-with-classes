#include "ValueNode.hpp"

ValueNode::ValueNode(std::size_t x_position
  , std::size_t y_position
  , std::size_t nx
  , double d_1
  , double d_2
  , bool b_1
  , int i_1)
  : Node(x_position, y_position, nx),
    v1 {{d_1, d_2}},
    b1 {b_1},
    i1 {i_1}
{}
