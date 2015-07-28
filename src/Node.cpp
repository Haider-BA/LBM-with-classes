#include "Node.hpp"

Node::Node(std::size_t x_position
  , std::size_t y_position
  , std::size_t nx)
  : x {x_position},
    y {y_position},
    n {y * nx + x},
    df_node {},
    neighbours {}
{}
