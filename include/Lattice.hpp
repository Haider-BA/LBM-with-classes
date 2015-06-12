#ifndef LATTICE_HPP_
#define LATTICE_HPP_
#include <vector>

class Lattice {
 public:
  /**
   * Constructor: (default) Override default constructor to set lattice to
   * zero
   */
  Lattice();

  /**
   * Constructor: Creates a 2D lattice with the specified dimensions [(num_rows
   * + 2) * (num_cols + 2)] to account for extra boundary layer. The depth of
   * each node is determined by the type of lattice (depth 9 for which_type == 0
   * and depth 2 for which_type == 1). Throws an exception for which_type values
   * other than 0 and 1.
   *
   * \param which_type specifies the type of lattice to be created
   * \param num_rows number of rows for the lattice
   * \param num_cols number of columns for the lattice
   * \throw 1 if which_type is not 0 or 1
   */
  Lattice(std::size_t num_dimensions
    , std::size_t num_discrete_velocities
    , std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt);

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
  std::size_t GetNumberOfDiscreteVelocities() const;

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

  /**
   * Initialize the 1D lattice with the same initial value
   * \param lattice the 1D lattice to be initialized
   * \param initial_value the value used to initialized the lattice
   */
  void Init(std::vector<double> &lattice
    , double initial_value);

  /**
   * Initialize the 2D lattice with the same initial values at each node
   * \param lattice the 2D lattice to be initialized
   * \param initial_values a 1D vector containing the initial values at each
   *        each node
   * \throw runtime_error if size of initial_values does not match depth of
   *        lattice
   */
  void Init(std::vector<std::vector<double>> &lattice
    , const std::vector<double> &initial_values);

  /**
   * Initialize the 2D lattice with the same values as another 2D lattice
   * \param lattice the 2D lattice to be initialized
   * \param initial_lattice the 2D lattice containing the initial lattice
   *        values
   * \throw runtime_error if size of initial_lattice does not match the size of
   *        lattice
   *        runtime_error if depth of initial_lattice does not match the depth
   *        of lattice
   */
  void Init(std::vector<std::vector<double>> &lattice
    , const std::vector<std::vector<double>> &initial_lattice);

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  void InitSrc(std::vector<double> &lattice_src
    , const std::vector<std::vector<unsigned>> &src_position
    , const std::vector<double> &src_strength);

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  void InitSrc(std::vector<std::vector<double>> &lattice_src
    , const std::vector<std::vector<int>> &src_position
    , const std::vector<std::vector<double>> &src_strength);

  /**
   * Calculates the equilibrium distribution function values at each node and
   * writes that to lattice_eq using <some formula>
   * \param lattice_eq the 2D lattice containing the equilibrium distribution
   *        function values
   * \param rho the density at each node (in lattice units)
   * \throw runtime_error if lattice_eq is not a depth 9 lattice
   */
  void ComputeEq(std::vector<std::vector<double>> &lattice_eq
    , const std::vector<double> &rho);

  /**
   * Flips the lattice for ease of printing out the lattice according to the
   * direction convention
   * \param lattice the lattice to be flipped
   * \return flipped lattice
   */
  std::vector<double> Flip(const std::vector<double> &lattice);

  /**
   * Flips the lattice for ease of printing out the lattice according to the
   * direction convention
   * \param lattice the lattice to be flipped
   * \return flipped lattice
   */
  std::vector<std::vector<double>> Flip(
      const std::vector<std::vector<double>> &lattice);

  /**
   * Prints the lattice for debugging purposes.
   * \param lattice the lattice to be printed
   */
  void Print(const std::vector<double> &lattice);

  /**
   * Prints the lattice for debugging purposes.
   * \param which_to_print selector to change way the lattice is printed
   *        0 each node is printed as a 3 * 3 square
   *        1 each node is printed as a single line
   */
  void Print(int which_to_print
    , const std::vector<std::vector<double>> &lattice);

  /**
   * NS distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f_;

  /**
   * CDE distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g_;

  /**
   * NS equilibrium distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f_eq_;

  /**
   * CDE equilibrium distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g_eq_;

  /**
   * Body force values for NS stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> src_f_;

  /**
   * Source term for CDE stored row-wise in a 1D vector.
   */
  std::vector<double> src_g_;

  /**
   * Density for NS stored row-wise in a 1D vector.
   */
  std::vector<double> rho_f_;

  /**
   * Density for CDE stored row-wise in a 1D vector.
   */
  std::vector<double> rho_g_;

  /**
   * Lattice velocity stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> u_;

 private:
  // 6  2  5  ^
  //  \ | /   |
  // 3--0--1  |
  //  / | \   |
  // 7  4  8  +------->
  enum DiscreteDirections {
    E = 1,
    N,
    W,
    S,
    NE,
    NW,
    SW,
    SE,
  };
  const std::vector<std::vector<double>> e_ = {{0, 0},
                                            {1, 0}, {0, 1}, {-1, 0}, {0, -1},
                                            {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  const std::vector<double> omega_ = {16. / 36.,
                                      4. / 36., 4. / 36., 4. / 36., 4. / 36.,
                                      1. / 36., 1. / 36., 1. / 36., 1. / 36.};
  std::size_t number_of_dimensions_;
  std::size_t number_of_discrete_velocities_;
  std::size_t number_of_rows_;
  std::size_t number_of_columns_;
  double space_step_;
  double time_step_;
  double c_;
  double input_parameter_check_value_ = space_step_ * time_step_ *
      number_of_dimensions_ * number_of_discrete_velocities_ * number_of_rows_ *
      number_of_columns_;
};
#endif  // LATTICE_HPP_
