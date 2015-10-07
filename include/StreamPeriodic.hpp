#ifndef STREAM_PERIODIC_HPP_
#define STREAM_PERIODIC_HPP_
#include <vector>
#include "LatticeModel.hpp"
#include "StreamModel.hpp"

class StreamPeriodic: public StreamModel {
 public:
  /**
   * Constructor: Creates a periodic streaming model for the D2Q9 lattice model
   * \param lm lattice model which contains information on the number of rows,
   *        columns, dimensions, discrete directons and lattice velocity
   */
  StreamPeriodic(LatticeModel &lm);

  /**
   * Destructor
   */
  ~StreamPeriodic() = default;

  /**
   * Performs the streaming step based on "Introduction to Lattice Boltzmann
   * Methods" and assumes there are periodic boundary conditions on all edges
   * and corners.
   * \param df lattice distribution functions stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> Stream(
      const std::vector<std::vector<double>> &df);
};

#endif  // STREAM_PERIODIC_HPP_
