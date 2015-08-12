#ifndef LATTICE_D2Q9_HPP_
#define LATTICE_D2Q9_HPP_
#include <vector>
#include "LatticeModel.hpp"

class LatticeD2Q9: public LatticeModel {
 public:
  /**
   * Constructor: Create lattice model for D2Q9 with the same velocity at each
   * node
   * \param num_rows number of rows
   * \param num_cols number of columns
   * \param dx space step
   * \param dt time step
   * \param initial_velocity initial velocity of the lattice
   */
  LatticeD2Q9(std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt
    , const std::vector<double> &initial_velocity);

  /**
   * Constructor: Create lattice model for D2Q9 with variable velocity at each
   * node
   * \param num_rows number of rows
   * \param num_cols number of columns
   * \param dx space step
   * \param dt time step
   * \param initial_velocity initial velocity of the lattice
   */
  LatticeD2Q9(std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt
    , const std::vector<std::vector<double>> &initial_velocity);

  /**
   * Destructor
   */
  virtual ~LatticeD2Q9() = default;
};
#endif  // LATTICE_D2Q9_HPP_
