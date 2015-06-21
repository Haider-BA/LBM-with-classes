#include "LatticeModel.hpp"
#include <stdexcept>  // runtime_error

std::size_t LatticeModel::GetNumberOfDimensions() const
{
  return number_of_dimensions_;
}

std::size_t LatticeModel::GetNumberOfDirections() const
{
  return number_of_directions_;
}

std::size_t LatticeModel::GetNumberOfRows() const
{
  return number_of_rows_;
}

std::size_t LatticeModel::GetNumberOfColumns() const
{
  return number_of_columns_;
}

double LatticeModel::GetSpaceStep() const
{
  return space_step_;
}

double LatticeModel::GetTimeStep() const
{
  return time_step_;
}

bool LatticeModel::CheckInput()
{
  return number_of_dimensions_ == 0 || number_of_directions_ == 0 ||
      number_of_rows_ == 0 || number_of_columns_ == 0;
}
