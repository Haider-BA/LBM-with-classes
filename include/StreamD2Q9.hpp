#ifndef STREAM_D2Q9_HPP_
#define STREAM_D2Q9_HPP_
#include <vector>
#include "LatticeModel.hpp"
#include "StreamModel.hpp"

class StreamD2Q9: public StreamModel {
 public:
  StreamD2Q9(LatticeModel &lm);

  ~StreamD2Q9() = default;

  std::vector<std::vector<double>> Stream(
      const std::vector<std::vector<double>> &df);

};

#endif // STREAM_D2Q9_HPP_

