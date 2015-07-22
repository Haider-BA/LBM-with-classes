#include "ValueNode.hpp"

ValueNode::ValueNode(std::size_t x_position
  , std::size_t y_position
  , std::size_t nx
  , double d_1
  , double d_2
  , bool b_1)
  : Node(x_position, y_position, nx),
    d1 {d_1},
    d2 {d_2},
    d3 {},
    b1 {b_1}
{}
