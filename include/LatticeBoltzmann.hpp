#ifndef LATTICEBOLTZMANN_HPP_
#define LATTICEBOLTZMANN_HPP_
#include "LatticeModel.hpp"

class LatticeBoltzmann {
 public:
  /**
   * Constructor: (default) Override default constructor to throw exception and
   * force user to initialize the lattice with values when calling the
   * constructor.
   * \throw std::runtime_error if called
   */
  LatticeBoltzmann();

  LatticeBoltzmann(const LatticeModel &lm);
};
#endif  // LATTICEBOLTZMANN_HPP_
