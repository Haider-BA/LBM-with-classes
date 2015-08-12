#ifndef STREAM_MODEL_CPP_
#define STREAM_MODEL_CPP_
#include <vector>
#include "LatticeModel.hpp"

class StreamModel {
 public:
  /**
   * Constructor: Base class for stream models
   * \param lm lattice model which contains information on the number of rows,
   *        columns, dimensions, discrete directions and lattice velocity
   */
  StreamModel(LatticeModel &lm);

  /**
   * Virtual destruction since we are deriving from this class
   */
  virtual ~StreamModel() = default;

  /**
   * Pure virtual function for the streaming function
   * \param df lattice distribution function stored row-wise in a 2D vector
   */
  virtual std::vector<std::vector<double>> Stream(
      const std::vector<std::vector<double>> &df) = 0;

 protected:
  /**
   * Enumeration for discrete directions to be used with distribution functions
   */
  enum Directions {
    E = 1,
    N,
    W,
    S,
    NE,
    NW,
    SW,
    SE
  };

  /**
   * Reference to lattice model
   */
  LatticeModel &lm_;
};
#endif // STREAM_MODEL_CPP_
