#ifndef NODE_HPP_
#define NODE_HPP_
#include <vector>
class Node {
 public:
  Node(std::size_t x_position
    , std::size_t y_position
    , std::size_t nx);
  virtual ~Node() = default;
  std::size_t x;
  std::size_t y;
  std::size_t n;
};
#endif // NODE_HPP_
