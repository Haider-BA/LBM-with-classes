#ifndef LATTICEMODEL_HPP_
#define LATTICEMODEL_HPP_
#include <vector>

class LatticeModel {
 public:
  /**
   * Constructor: creates lattice model
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
   * Constructor: creates lattice model
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
   // temporary workaround due to MinGW bug
  virtual ~LatticeModel(){};// = default;

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
   * Pure virtual function that compute the density of the lattice, derived
   * classes will use their own get moments function to compute density
   * \param lattice 2D vector containing distribution functions
   * \return density of the lattice stored row-wise
   */
  virtual std::vector<double> ComputeRho(
      const std::vector<std::vector<double>> &lattice) = 0;

  /**
   * Pure virtual function that compute the velocity of the lattice, derived
   * classes will use their own get moments function to compute velocity
   * \param lattice 2D vector containing distribution functions
   * \param rho density of the lattice stored row-wise
   * \param src source term of the lattice stored row-wise
   * \return density of the lattice stored row-wise
   */
  virtual std::vector<std::vector<double>> ComputeU(
      const std::vector<std::vector<double>> &lattice
    , const std::vector<double> &rho
    , const std::vector<std::vector<double>> &src) = 0;

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
  bool CheckInput();
  std::size_t number_of_dimensions_;
  std::size_t number_of_directions_;
  std::size_t number_of_rows_;
  std::size_t number_of_columns_;
  double space_step_;
  double time_step_;
  double c_ = space_step_ / time_step_;
};
#endif  // LATTICEMODEL_HPP_
