#ifndef LATTICEMODEL_HPP_
#define LATTICEMODEL_HPP_
#include <vector>

class LatticeModel {
 public:
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
   * Get the number of rows of the lattice, will be 2 smaller than the
   * actual number of rows created due to additional boundary layer.
   * \return number of rows of the lattice
   */
  std::size_t GetNumberOfRows() const;

  /**
   * Get the number of columns of the lattice, will be 2 smaller than
   * the actual number of columns created due to additional boundary layer.
   * \return number of columns of the lattice
   */
  std::size_t GetNumberOfColumns() const;
  double GetSpaceStep() const;
  double GetTimeStep() const;

 protected:
  bool CheckInput();
  std::size_t number_of_dimensions_;
  std::size_t number_of_directions_;
  std::size_t number_of_rows_;
  std::size_t number_of_columns_;
  double space_step_;
  double time_step_;
};
#endif  // LATTICEMODEL_HPP_
