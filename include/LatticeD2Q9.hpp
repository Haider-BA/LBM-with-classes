#ifndef LATTICED2Q9_HPP_
#define LATTICED2Q9_HPP_
#include "LatticeModel.hpp"
#include <vector>

class LatticeD2Q9: public LatticeModel {
 public:
  LatticeD2Q9(std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt);

  ~LatticeD2Q9() = default;

  /**
   * Calculates zeroth moment of a node based on formula in LBIntro
   * \param node distribution function node containing nine discrete velocity
   *        vectors
   * \return zeroth moment of node
   */
  double GetZerothMoment(const std::vector<double> &node);

  /**
   * Calculates first moment of a node based on formula in LBIntro
   * \param node distribution function node containing nine discrete velocity
   *        vectors
   * \return first moment of node
   */
  std::vector<double> GetFirstMoment(const std::vector<double> &node);

  /** \brief
   *
   * \param
   * \param
   * \return
   *
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
};
#endif  // LATTICED2Q9_HPP_
