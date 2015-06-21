#ifndef LATTICEBOLTZMANN_HPP_
#define LATTICEBOLTZMANN_HPP_
#include <vector>
#include "LatticeModel.hpp"

class LatticeBoltzmann {
 public:
  /**
   * Constructor: (default) Override default constructor to throw exception to
   * to forbid declaring uninitialized
   * Does not work because of uninitialized reference member error (lm_)
   */
//  LatticeBoltzmann();

  /**
   * Constructor: Creates lattice
   * \param
   * \param
   * \return
   *
   */
  LatticeBoltzmann(double t_total
    , double diffusion_coefficient
    , double kinematic_viscosity
    , double initial_density_f
    , double initial_density_g
    , const std::vector<double> &u0
    , const std::vector<std::vector<std::size_t>> &src_position_f
    , const std::vector<std::vector<double>> &src_strength_f
    , const std::vector<std::vector<std::size_t>> &src_position_g
    , const std::vector<double> &src_strength_g
    , const std::vector<std::vector<std::size_t>> &obstacles_pos
    , bool is_ns
    , bool is_cd
    , bool is_instant
    , bool has_obstacles
    , LatticeModel &lm);

   /**
   * Get the number of dimensions of the lattice. 2 for 2D and 3 for 3D.
   *
   * \return number of dimensions of the lattice
   */
  std::size_t GetNumberOfDimensions() const;

  /**
   * Get the number of discrete velocities of the lattice, specified by the
   * model used. 9 for Q9.
   *
   * \return number of discrete velocities of the lattice
   */
  std::size_t GetNumberOfDirections() const;

  /**
   * Get the number of rows of the lattice, will be 2 smaller than the
   * actual number of rows created due to additional boundary layer.
   *
   * \return number of rows of the lattice
   */
  std::size_t GetNumberOfRows() const;

  /**
   * Get the number of columns of the lattice, will be 2 smaller than
   * the actual number of columns created due to additional boundary layer.
   *
   * \return number of columns of the lattice
   */
  std::size_t GetNumberOfColumns() const;

  void Init(std::vector<bool> &lattice
  , const std::vector<std::vector<std::size_t>> &position);

  std::vector<double> Flip(const std::vector<double> &lattice);
  std::vector<std::vector<double>> Flip(
    const std::vector<std::vector<double>> &lattice);
  void Print(const std::vector<double> &lattice);
  void Print(int which_to_print
  , const std::vector<std::vector<double>> &lattice);

  /**
   * Lattice velocity stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> u;

  /**
   * NS distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f;

  /**
   * Density for NS stored row-wise in a 1D vector.
   */
  std::vector<double> rho_f;

  /**
   * Boundary nodes for NS lattice
   */
  std::vector<std::vector<double>> boundary_f;

  /**
   * CDE distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g;

  /**
   * Density for CDE stored row-wise in a 1D vector.
   */
  std::vector<double> rho_g;

  /**
   * Boundary nodes for CDE lattice
   */
  std::vector<std::vector<double>> boundary_g;


  /**
   * Lattice containing obstacles
   */
  std::vector<bool> obstacles;

 private:
  /**
   * Checks input parameters to ensure there's not invalid values
   * \return true if there is invalid: 0 in any values, both is_ns_ and is_cd_
   *         are false
   *         false if all input values are valid
   */
  bool CheckParameters();
  // input parameters
  double total_time_;
  bool is_ns_;
  bool is_cd_;
  bool is_instant_;
  bool has_obstacles_;
  // LatticeModel to take care of dims, dirs, rows, cols and discrete e vectors
  // by reference, similar to by pointer
  // https://stackoverflow.com/questions/9285627/is-it-possible-to-pass-derived-
  // classes-by-reference-to-a-function-taking-base-cl
  LatticeModel &lm_;
};
#endif // LATTICEBOLTZMANN_HPP_
