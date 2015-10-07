#ifndef BOUNCE_BACK_NODES_HPP_
#define BOUNCE_BACK_NODES_HPP_
#include <vector>
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "Node.hpp"
#include "StreamModel.hpp"

class BouncebackNodes: public BoundaryNodes {
 public:
  /**
   * Creates a full-way bounceback nodes according to "http://lbmworkshop.com/wp
   * -content/uploads/2011/08/Straight_boundaries.pdf"
   * \param lm LatticeModel to provide information on number of rows, columns,
   *        dimensions, discrete directions and lattice velocity
   * \param cm CollisionModel to indicate which nodes to be skipped during the
   *        collision step
   */
  BouncebackNodes(LatticeModel &lm
    , CollisionModel *cm);

  /**
   * Creates a half-way bounceback node according to "http://lbmworkshop.com/wp-
   * content/uploads/2011/08/Straight_boundaries.pdf"
   * \param lm LatticeModel to provide information on number of rows, columns,
   *        dimensions, discrete directions and lattice velocity
   * \param sm StreamModel to copy the prestream node distribution functions
   *        so they can be bounced back in the same time step
   */
  BouncebackNodes(LatticeModel &lm
    , StreamModel *sm);

  /**
   * Override copy constructor due to -Weffc++ warnings
   */
  BouncebackNodes(const BouncebackNodes&) = default;

  /**
   * Override copy assignment due to -Weffc++ warnings
   */
  BouncebackNodes& operator= (const BouncebackNodes&) = default;

  /**
   * Destructor
   */
  ~BouncebackNodes() = default;

  /**
   * Adds a bounceback node
   * \param x x-coordinate of the node
   * \param y y-coordinate of the node
   */
  void AddNode(std::size_t x, std::size_t y);

  /**
   * Performs the bounceback boundary condition on the boundary nodes based on
   * the type of bounceback nodes used.
   * Full-way bounceback: Reflects all the node distribution functions in the
   *     opposite direction (except center distribution function)
   * Half-way bounceback: Copies the prestream node distribution functions
   *    before streaming. Updates the post-stream unknown distribution functions
   *    with the prestream distribution functions in the opposite directions
   * \param df lattice distribution functions store row-wise in a 2D vector
   * \param is_modify_stream Boolean toggle for half-way bounceback as it has
   *        both pre-stream and post-stream functions. Used to fit in with how
   *        all boundary conditions are called
   */
  void UpdateNodes(std::vector<std::vector<double>> &df
    , bool is_modify_stream);

  /**
   * Vector used to store information about the boundary nodes such as their
   * position in the lattice
   */
  std::vector<Node> nodes;

 protected:
  /**
   * Pointer to collision model as half-way bounceback nodes do not require
   * collision models and NULL references can't be declared
   */
  CollisionModel* cm_ = nullptr;

  /**
   * Pointer to stream model as full-way bounceback nodes do not require stream
   * models and NULL references can't be declared
   */
  StreamModel* sm_ = nullptr;
};
#endif  // BOUNCE_BACK_NODES_HPP_
