#ifndef LATTICEBOLTZMANN_HPP_
#define LATTICEBOLTZMANN_HPP_
#include <vector>

class LatticeBoltzmann {
 public:
  /**
   * Constructor: (default) Override default constructor to set lattice to
   * zero
   */
  LatticeBoltzmann();

  LatticeBoltzmann(std::size_t num_dims
    , std::size_t num_dirs
    , std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt
    , double t_total
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
    , bool has_obstacles);

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

  /**
   * Lattice velocity stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> u;

  /**
   * NS distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f;

  /**
   * NS equilibrium distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f_eq;

  /**
   * Density for NS stored row-wise in a 1D vector.
   */
  std::vector<double> rho_f;

  /**
   * Boundary nodes for NS lattice
   */
  std::vector<std::vector<double>> boundary_f;

  /**
   * Body force values for NS stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> src_f;

  /**
   * CDE distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g;

  /**
   * CDE equilibrium distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g_eq;

  /**
   * Density for CDE stored row-wise in a 1D vector.
   */
  std::vector<double> rho_g;

  /**
   * Boundary nodes for CDE lattice
   */
  std::vector<std::vector<double>> boundary_g;

  /**
   * Source term for CDE stored row-wise in a 1D vector.
   */
  std::vector<double> src_g;

  /**
   * Lattice containing obstacles
   */
  std::vector<bool> obstacles;

 private:
  bool CheckParameters();
  // input parameters
  std::size_t number_of_dimensions_;
  std::size_t number_of_directions_;
  std::size_t number_of_rows_;
  std::size_t number_of_columns_;
  double dx_;
  double dt_;
  double total_time_;
  bool is_ns_;
  bool is_cd_;
  bool is_instant_;
  bool has_obstacles_;
  // calculated parameters
  double c_;
  double cs_sqr_;
  double tau_ns_;
  double tau_cd_;
};
#endif // LATTICEBOLTZMANN_HPP_
