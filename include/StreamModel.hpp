#ifndef STREAM_MODEL_CPP_
#define STREAM_MODEL_CPP_
#include <vector>
#include "LatticeModel.hpp"

class StreamModel {
 public:
  StreamModel(LatticeModel &lm);

  virtual ~StreamModel() = default;

  void AddNodeToBounceback(std::size_t n);

  virtual std::vector<std::vector<double>> Stream(
      const std::vector<std::vector<double>> &df) = 0;

  std::vector<bool> bounce_back;

 protected:
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

  LatticeModel &lm_;
};
#endif // STREAM_MODEL_CPP_
