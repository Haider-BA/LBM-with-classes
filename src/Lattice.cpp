#include "Lattice.hpp"
#include <vector>
#include <iostream>
#include <iomanip> // std::setprecision

Lattice::Lattice()
  : number_of_rows_ {0},
    number_of_columns_ {0}
{}

Lattice::Lattice(int which_type, std::size_t num_rows, std::size_t num_cols)
  : lattice_type_ {which_type},
    number_of_rows_ {num_rows},
    number_of_columns_ {num_cols}
{
  // 0 for distribution functions, 1 for velocities,
  switch (which_type) {
    case 0:
      values_.assign((num_rows + 2) * (num_cols + 2),
                     std::vector<double>(number_of_discrete_velocities_));
      break;
    case 1:
      values_.assign((num_rows + 2) * (num_cols + 2),
                     std::vector<double>(number_of_dimensions_));
      break;
    default:
      throw 1;
  }
}

std::size_t Lattice::GetNumberOfDimensions() const
{
  return number_of_dimensions_;
}

std::size_t Lattice::GetNumberOfDiscreteVelocities() const
{
  return number_of_discrete_velocities_;
}

std::size_t Lattice::GetNumberOfRows() const
{
  return number_of_rows_;
}

std::size_t Lattice::GetNumberOfColumns() const
{
  return number_of_columns_;
}

void Lattice::Print()
{
  switch (lattice_type_) {
    case 0:
      // row of lattice
      for (auto y = 0; y < (number_of_columns_ + 2) * (number_of_rows_ + 2);
           y += number_of_columns_ + 2) {
        // rows in the Q9 square
        for (auto i = 0; i < 3; ++i) {
          // column of lattice
          for (auto x = 0; x < number_of_columns_ + 2; ++x) {
            if (i == 0){
            std::cout << std::fixed << std::setprecision(2)
                      << values_[y + x][6] << " "
                      << values_[y + x][2] << " "
                      << values_[y + x][5] << " ";
            }
            else if (i == 1){
            std::cout << std::fixed << std::setprecision(2)
                      << values_[y + x][3] << " "
                      << values_[y + x][0] << " "
                      << values_[y + x][1] << " ";
            }
            else if (i == 2){
            std::cout << std::fixed << std::setprecision(2)
                      << values_[y + x][7] << " "
                      << values_[y + x][4] << " "
                      << values_[y + x][8] << " ";
            }
            std::cout << "  ";
          } // x
          std::cout << "\n";
        } // i
        std::cout << "\n";
      } // y
      break;
    case 1:
      // row of lattice
      for (auto y = 0; y < (number_of_columns_ + 2) * (number_of_rows_ + 2);
           y += number_of_columns_ + 2) {
        // column of lattice
        for (auto x = 0; x < number_of_columns_ + 2; ++x) {
          // elements of one line of Q9 square
          for(auto d = 0; d < 2; ++d) {
            std::cout << values_[y + x][d] << " ";
          } // d
          std::cout << "  ";
        } // x
        std::cout << "\n";
      } // y
      break;
    default:
      std::cout << "Not a 2D lattice\n";
      break;
  }
}
