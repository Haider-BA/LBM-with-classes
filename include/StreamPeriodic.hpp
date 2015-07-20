#ifndef STREAM_PERIODIC_HPP_
#define STREAM_PERIODIC_HPP_
#include <vector>
#include "LatticeModel.hpp"
#include "StreamModel.hpp"

class StreamPeriodic: public StreamModel {
 public:
  StreamPeriodic(LatticeModel &lm);

  ~StreamPeriodic() = default;

  std::vector<std::vector<double>> Stream(
      const std::vector<std::vector<double>> &df);
};

#endif // STREAM_PERIODIC_HPP_
