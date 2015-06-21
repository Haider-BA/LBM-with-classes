#ifndef LATTICED2Q9_HPP_
#define LATTICED2Q9_HPP_
#include "LatticeModel.hpp"

class LatticeD2Q9: public LatticeModel {
 public:
  LatticeD2Q9(std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt);

  std::vector<std::vector<double>> c_ = {{0, 0},
                                         {1, 0}, {0, 1}, {-1, 0}, {0, -1},
                                         {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  std::vector<double> omega_ = {16.0 / 36.0,
                                4.0 / 36.0, 4.0 / 36.0, 4.0 / 36.0, 4.0 / 36.0,
                                1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
};
#endif  // LATTICED2Q9_HPP_
