#ifndef NODE_HPP_
#define NODE_HPP_
#include <vector>
class Node {
 public:
  /**
   * Constructor: Creates a node which stores information about its position in
   * the lattice
   * \param x_position x-coordinate of the node
   * \param y_position y-coordinate of the node
   * \param nx number of columns of the lattice, for calculating index in the
   *        distribution function vector
   */
  Node(std::size_t x_position
    , std::size_t y_position
    , std::size_t nx);

  /**
   * Virtual destructor since we are deriving from this class
   */
  virtual ~Node() = default;

  /**
   * x-coordinate
   */
  std::size_t x;

  /**
   * y-coordinate
   */
  std::size_t y;

  /**
   * Index in distribution function vector
   */
  std::size_t n;

  /**
   * Distribution functions of a single node, used by half-way bounceback nodes
   */
  std::vector<double> df_node;
};
#endif  // NODE_HPP_
