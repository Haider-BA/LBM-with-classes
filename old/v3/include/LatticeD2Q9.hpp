#ifndef LATTICED2Q9_HPP_
#define LATTICED2Q9_HPP_
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

  /**
   * Compute density at each node by summing up its distribution functions
   * \param lattice 2D vector containing distribution functions
   * \return density of lattice stored row-wise in a 1D vector
   */
  std::vector<double> ComputeRho(
      const std::vector<std::vector<double>> &lattice);

  /**
   * Calculated velocity for NS equation based on formula in Guo2002
   * \param lattice 2D vector containing distribution functions of the NS
   *        equation
   * \param rho 1D vector containing the density at each node of lattice
   * \param src 2D vector containing body force density at each node of lattice
   * \return 2D vector containing velocity at each node of lattice
   */
  std::vector<std::vector<double>> ComputeU(
      const std::vector<std::vector<double>> &lattice
    , const std::vector<double> &rho
    , const std::vector<std::vector<double>> &src);

  std::vector<std::vector<double>> ComputeULid(
      const std::vector<std::vector<double>> &lattice
    , const std::vector<double> &rho
    , const std::vector<std::vector<double>> &src
    , double u_lid);

//  std::vector<std::vector<double>> e_d2q9 = {{0, 0},
//      {1, 0}, {0, 1}, {-1, 0}, {0, -1},
//      {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
};
#endif  // LATTICED2Q9_HPP_
