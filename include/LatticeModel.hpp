#ifndef LATTICEMODEL_HPP_
#define LATTICEMODEL_HPP_
#include <vector>

class LatticeModel {
 public:
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
  double c_;
};
#endif  // LATTICEMODEL_HPP_
