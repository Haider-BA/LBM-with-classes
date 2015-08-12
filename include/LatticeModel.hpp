#ifndef LATTICE_MODEL_HPP_
#define LATTICE_MODEL_HPP_
#include <vector>

class LatticeModel {
 public:
  /**
   * Constructor: creates lattice model with the same velocity at each node
   * \param num_dims number of dimensions
   * \param num_dirs number of discrete directions
   * \param num_rows number of rows
   * \param num_cols number of columns
   * \param dx space step
   * \param dt time step
   * \param initial_velocity initial velocity of the lattice
   */
  LatticeModel(std::size_t num_dims
    , std::size_t num_dirs
    , std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt
    , const std::vector<double> &initial_velocity);

  /**
   * Constructor: creates lattice model with variable velocity at each node
   * \param num_dims number of dimensions
   * \param num_dirs number of discrete directions
   * \param num_rows number of rows
   * \param num_cols number of columns
   * \param dx space step
   * \param dt time step
   * \param initial_velocity initial velocity of the lattice
   */
  LatticeModel(std::size_t num_dims
    , std::size_t num_dirs
    , std::size_t num_rows
    , std::size_t num_cols
    , double dx
    , double dt
    , const std::vector<std::vector<double>> &initial_velocity);

  /**
   * Virtual destructor since we deriving from this class, see Collision.hpp
   */
  virtual ~LatticeModel() = default;

  /**
   * Get the number of dimensions of the lattice. 2 for 2D and 3 for 3D.
   * \return number of dimensions of the lattice
   */
  std::size_t GetNumberOfDimensions() const;

  /**
   * Get the number of discrete velocities of the lattice, specified by the
   * model used. 9 for Q9.
   * \return number of discrete velocities of the lattice
   */
  std::size_t GetNumberOfDirections() const;

  /**
   * Get the number of rows of the lattice
   * \return number of rows of the lattice
   */
  std::size_t GetNumberOfRows() const;

  /**
   * Get the number of columns of the lattice
   * \return number of columns of the lattice
   */
  std::size_t GetNumberOfColumns() const;

  /**
   * Get the space step (dx) of the model
   * \return space step of the model
   */
  double GetSpaceStep() const;

  /**
   * Get the time step (dt) of the model
   * \return time step of the model
   */
  double GetTimeStep() const;

  /**
   * Get the lattice speed (c) of the model
   * \return lattice speed of the model
   */
  double GetLatticeSpeed() const;

  /**
   * Lattice velocity stored row-wise in a 2D vector
   */
  std::vector<std::vector<double>> u;

  /**
   * Unit vector for each discrete velocity direction
   */
  std::vector<std::vector<double>> e;

  /**
   * Weights each discrete velocity according to LBIntro
   */
  std::vector<double> omega;

 protected:
  /**
   * Checks if input parameters for lattice model is valid, prevents creation of
   * invalid lattice model, such as a size 0 x 0 lattice
   * \return validity of lattice model
   *         TRUE: input parameters are invalid
   *         FALSE: input parameters are valid
   */
  bool CheckInput();

  /**
   * Number of dimensions of the lattice model, 2 for D2Q9 model
   */
  std::size_t number_of_dimensions_;

  /**
   * Number of discrete directions of the lattice model, 9 for D2Q9 model
   */
  std::size_t number_of_directions_;

  /**
   * Number of rows of the lattice model
   */
  std::size_t number_of_rows_;

  /**
   * Number of columns of the lattice model
   */
  std::size_t number_of_columns_;

  /**
   * Space step of the lattice model, dx
   */
  double space_step_;

  /**
   * Time step of the lattice model, dt
   */
  double time_step_;

  /**
   * Propagation speed on the lattice. Based on "Introduction to Lattice
   * Boltzmann Methods"
   */
  double c_ = space_step_ / time_step_;
};
#endif  // LATTICE_MODEL_HPP_
