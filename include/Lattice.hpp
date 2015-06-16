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
   * \param dx
   * \param dt
   * \param diffusion_coefficient
   * \param kinematic_viscosity
   * \param is_cd
   * \param is_ns
   * \param is_instant
   * \throw std::runtime_error if there is incorrect values in the input
   *        parameters
   */
  Lattice(std::size_t num_dimensions
    , std::size_t num_discrete_velocities
    , std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt
    , double t_total
    , double diffusion_coefficient
    , double kinematic_viscosity
    , double density_f
    , double density_g
    , const std::vector<double> &u0
    , const std::vector<std::vector<unsigned>> &src_pos_f
    , const std::vector<std::vector<unsigned>> &src_pos_g
    , const std::vector<std::vector<double>> &src_strength_f
    , const std::vector<double> &src_strength_g
    , bool is_cd
    , bool is_ns
    , bool is_instant);

  Lattice(const Lattice&) = default;

  ~Lattice() = default;

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

  /**
   * Initialize a depth 1 2D source lattice by reading from the source position
   * vectors and the source magnitude vectors,
   * \param lattice_src 2D lattice containing the source strength of each node
   * \param src_position 2D vector containing position information about the
   *        source
   * \param src_strength 2D vector containing strength information about the
   *        source
   * \throw runtime_error if size of src_position does not match size of
   *        src_strength
   * \throw runtime_error if size of vector in src_position does not match depth
   *        of lattice
   * \throw runtime_error if component in vector in src_position exceeds number
   *        of columns or rows of lattice
   */
  void InitSrc(std::vector<double> &lattice_src
    , const std::vector<std::vector<unsigned>> &src_position
    , const std::vector<double> &src_strength);

  /**
   * Initialize a depth 2 2D source lattice by reading from the source position
   * vectors and the source magnitude vectors,
   * \param lattice_src 2D lattice containing the source strength of each node
   * \param src_position 2D vector containing position information about the
   *        source
   * \param src_strength 2D vector containing strength information about the
   *        source
   * \throw runtime_error if size of src_position does not match size of
   *        src_strength
   * \throw runtime_error if size of vector in src_position does not match depth
   *        of lattice
   * \throw runtime_error if component in vector in src_position exceeds number
   *        of columns or rows of lattice
   */
  void InitSrc(std::vector<std::vector<double>> &lattice_src
    , const std::vector<std::vector<unsigned>> &src_position
    , const std::vector<std::vector<double>> &src_strength);

  /** \brief
   *
   * \return void
   *
   */
  void InitAll();

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
   * 2D vector to store distribution functions values for nodes on the edges of
   * the lattice to be used in the streaming step. Current implementations are:
   * periodic boundary for the sides, no-slip boundary for top and bottom,
   * bounce-back for corners
   * \param lattice 2D lattice containing distribution functions at each node
   */
  void BoundaryCondition(const std::vector<std::vector<double>> &lattice
    , std::vector<std::vector<double>> &boundary);

  /**
   * Streams the lattice based on LBIntro. Takes node values from the boundary
   * 2D vector for nodes at the edges of lattice.
   * Chirila, D. B., (2010). Introduction to Lattice Boltzmann Methods
   * \param lattice 2D lattice containing pre-stream values of distribution
   *        functions
   * \param boundary 2D vector containing pre-stream values for nodes on the
   *        edges of the lattice badges on the pre-defined boundary conditions
   * \return 2D lattice containing post-stream values of distribution functions
   */
  std::vector<std::vector<double>> Stream(
      const std::vector<std::vector<double>> &lattice
    , const std::vector<std::vector<double>> &boundary);

  /**
   * Performs the collision step based on "A new scheme for source term in LBGK
   * model for convection–diffusion equation"
   * \param lattice 2D lattice containing the distribution functions
   * \param lattice_eq 2D lattice containing the equilibrium distribution
   *        functions
   */
  void Collide(std::vector<std::vector<double>> &lattice
    , const std::vector<std::vector<double>> &lattice_eq
    , std::vector<double> &src);

  /**
   * Performs the collision step based on Guo2002
   * \param lattice 2D lattice containing the distribution functions
   * \param lattice_eq 2D lattice containing the equilibrium distribution
   *        functions
   * \param src 2D lattice containing the source magnitude at each node
   * \param rho 2D lattice containing the density at each node
   *
   */
  void Collide(std::vector<std::vector<double>> &lattice
    , const std::vector<std::vector<double>> &lattice_eq
    , std::vector<std::vector<double>> &src
    , const std::vector<double> &rho);

  /**
   * Computes density at each node of the lattice
   * \param lattice 2D lattice containing the distribution functions
   * \return 2D lattice containing the density value at each node
   */
  std::vector<double> ComputeRho(
      const std::vector<std::vector<double>> &lattice);

  /**
   * Computes macroscopic velocity of the lattice
   * \param lattice 2D lattice containing the distribution functions of the NS
   *        equation
   * \param rho lattice containing the density of the NS equation
   * \param src lattice containing the source of the NS equation
   * \throw std::runtime_error if the size of the 3 input lattices do not match
   */
  void ComputeU(const std::vector<std::vector<double>> &lattice
  , const std::vector<double> &rho
  , const std::vector<std::vector<double>> &src);

  /**
   * Performance a series of streaming and collision depending on is_cd and
   * is_ns toggles
   */
  void TakeStep();

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
   *        2 for boundary
   *        3 for lattice with boundary
   */
  void Print(int which_to_print
    , const std::vector<std::vector<double>> &lattice);

  /**
   * NS distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f;

  /**
   * CDE distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g;

  /**
   * NS equilibrium distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f_eq;

  /**
   * CDE equilibrium distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g_eq;

  /**
   * Body force values for NS stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> src_f;

  /**
   * Source term for CDE stored row-wise in a 1D vector.
   */
  std::vector<double> src_g;

  /**
   * Density for NS stored row-wise in a 1D vector.
   */
  std::vector<double> rho_f;

  /**
   * Density for CDE stored row-wise in a 1D vector.
   */
  std::vector<double> rho_g;

  /**
   * Boundary nodes for NS lattice
   */
  std::vector<std::vector<double>> boundary_f;

  /**
   * Boundary nodes for CDE lattice
   */
  std::vector<std::vector<double>> boundary_g;

  /**
   * Lattice velocity stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> u;

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
    SE
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
  double total_time_;
  double diffusion_coefficient_;
  double kinematic_viscosity_;
  double c_;
  double cs_sqr_;
  double tau_cd_;
  double tau_ns_;
  double initial_density_f_;
  double initial_density_g_;
  std::vector<double> initial_velocity_;
  std::vector<std::vector<unsigned>> source_position_f_;
  std::vector<std::vector<unsigned>> source_position_g_;
  std::vector<std::vector<double>> source_strength_f_;
  std::vector<double> source_strength_g_;
  bool is_cd_;
  bool is_ns_;
  bool is_instant_;
  double input_parameter_check_value_ = space_step_ * time_step_ * total_time_ *
      number_of_dimensions_ * number_of_discrete_velocities_ * number_of_rows_ *
      initial_density_f_ * initial_density_g_ * initial_velocity_.size() *
      number_of_columns_ * diffusion_coefficient_ * kinematic_viscosity_ *
      (is_cd_ || is_ns_);
};
#endif  // LATTICE_HPP_
