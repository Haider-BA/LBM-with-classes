#ifndef LATTICE_HPP_
#define LATTICE_HPP_
#include <vector>

class Lattice {
 public:
  /**
   * Constructor: (default) Override default constructor to set lattice to
   * zero
   *
   */
  Lattice();

  /**
   * Constructor: Creates a 2D lattice with the specified dimensions [(num_rows
   * + 2) * (num_cols + 2)] to account for extra boundary layer. The depth of
   * each node is determined by the type of lattice (depth 9 for which_type == 0
   * and depth 2 for which_type == 1). Throws an exception for which_type values
   * other than 0 and 1.
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

  /** \brief Get the number of dimensions of the lattice. 2 for 2D and 3 for 3D.
   *
   * \return number of dimensions of the lattice
   *
   */
  std::size_t GetNumberOfDimensions() const;

  /** \brief Get the number of discrete velocities of the lattice, specified by
   * the model used. 9 for Q9.
   *
   * \return number of discrete velocities of the lattice
   *
   */
  std::size_t GetNumberOfDiscreteVelocities() const;

  /** \brief Get the number of rows of the lattice, will be 2 smaller than the
   * actual number of rows created due to additional boundary layer.
   *
   * \return number of rows of the lattice
   *
   */
  std::size_t GetNumberOfRows() const;

  /** \brief Get the number of columns of the lattice, will be 2 smaller than
   * the actual number of columns created due to additional boundary layer.
   *
   * \return number of columns of the lattice
   *
   */
  std::size_t GetNumberOfColumns() const;

  /** \brief
   *
   * \param initial_u const std::vector<double>
   * \return void
   *
   */
  void Init(std::vector<double> &lattice
    , double initial_value);

  /** \brief
   *
   * \param initial_u const std::vector<double>
   * \return void
   *
   */
  void Init(std::vector<std::vector<double>> &lattice
    , const std::vector<double> &initial_values);

  /** \brief
   *
   * \param initial_u const std::vector<double>
   * \return void
   *
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
  void ComputeEq(std::vector<std::vector<double>> &lattice_eq
    , const std::vector<double> &rho);

  /** \brief Prints the lattice for debugging purposes. Does not require
   * additional parameters as it will print based on the value of lattice_type_
   *
   * \return void
   *
   */
  void Print(const std::vector<double> &lattice);

  /** \brief Prints the lattice for debugging purposes. Does not require
   * additional parameters as it will print based on the value of lattice_type_
   *
   * \return void
   *
   */
  void Print(int which_to_print
    , const std::vector<std::vector<double>> &lattice);

  /**
   * Lattice entries stored row-wise in a 2D vector.
   *
   */
  std::vector<std::vector<double>> f_;
  std::vector<std::vector<double>> g_;
  std::vector<std::vector<double>> f_eq_;
  std::vector<std::vector<double>> g_eq_;
  std::vector<std::vector<double>> src_f_;
  std::vector<double> src_g_;
  std::vector<double> rho_f_;
  std::vector<double> rho_g_;
  std::vector<std::vector<double>> u_;

 private:
  // 6  2  5
  //  \ | /
  // 3--0--1
  //  / | \
  // 7  4  8
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
