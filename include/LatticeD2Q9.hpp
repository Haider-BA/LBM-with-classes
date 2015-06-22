#ifndef LATTICED2Q9_HPP_
#define LATTICED2Q9_HPP_
#include "LatticeModel.hpp"

class LatticeD2Q9: public LatticeModel {
 public:
  LatticeD2Q9(std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt);

  ~LatticeD2Q9() = default;
};
#endif  // LATTICED2Q9_HPP_
