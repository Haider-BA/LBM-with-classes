#ifndef STREAM_D2Q9_HPP_
#define STREAM_D2Q9_HPP_
#include <vector>
#include "LatticeModel.hpp"
#include "StreamModel.hpp"

class StreamD2Q9: public StreamModel {
 public:
  /**
   * Constructor: Creates a non-periodic streaming model for D2Q9 lattice model
   * \param lm Lattice model which contains information on the number of rows,
   *        columns, dimensions, discrete directions and lattice velocity
   */
  StreamD2Q9(LatticeModel &lm);

  /**
   * Destructor
   */
  ~StreamD2Q9() = default;

  /**
   * Performs the streaming function based on "Introduction to Lattice Boltzmann
   * Methods". Distribution functions which require off-lattice streaming are
   * unchanged
   * \param df lattice distribution functions stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> Stream(
      const std::vector<std::vector<double>> &df);
};

#endif // STREAM_D2Q9_HPP_

