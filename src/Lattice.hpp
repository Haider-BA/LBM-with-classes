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

  Lattice(int which_type, std::size_t num_rows, std::size_t num_cols);

  std::size_t GetNumberOfDimensions() const;

  std::size_t GetNumberOfDiscreteVelocities() const;

  std::size_t GetNumberOfRows() const;

  std::size_t GetNumberOfColumns() const;



  void Print();

  int lattice_type_ = -1;

  std::vector<std::vector<double>> values_;

 private:
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
  std::size_t number_of_dimensions_ = 2;
  std::size_t number_of_discrete_velocities_ = 9;
  std::size_t number_of_rows_;
  std::size_t number_of_columns_;
  const std::vector<std::vector<double>> e = {{0, 0},
                                            {1, 0}, {0, 1}, {-1, 0}, {0, -1},
                                            {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  const std::vector<double> omega = {16. / 36.,
                                     4. / 36., 4. / 36., 4. / 36., 4. / 36.,
                                     1. / 36., 1. / 36., 1. / 36., 1. / 36.};
};
#endif // LATTICE_HPP_
