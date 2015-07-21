#include "LatticeModel.hpp"
#include <stdexcept>  // runtime_error
#include <vector>

LatticeModel::LatticeModel(std::size_t num_dims
  , std::size_t num_dirs
  , std::size_t num_rows
  , std::size_t num_cols
  , double dx
  , double dt
  , const std::vector<double> &initial_velocity)
  : u {},
    e {},  // cannot pass in LatticeD2Q9 public member e_d2q9
    omega {},
    number_of_dimensions_ {num_dims},
    number_of_directions_ {num_dirs},
    number_of_rows_ {num_rows},
    number_of_columns_ {num_cols},
    space_step_ {dx},
    time_step_ {dt}
{
  u.assign(num_rows * num_cols, initial_velocity);
}

LatticeModel::LatticeModel(std::size_t num_dims
  , std::size_t num_dirs
  , std::size_t num_rows
  , std::size_t num_cols
  , double dx
  , double dt
  , const std::vector<std::vector<double>> &initial_velocity)
  : u {initial_velocity},
    e {},
    omega {},
    number_of_dimensions_ {num_dims},
    number_of_directions_ {num_dirs},
    number_of_rows_ {num_rows},
    number_of_columns_ {num_cols},
    space_step_ {dx},
    time_step_ {dt}
{}

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

double LatticeModel::GetLatticeSpeed() const
{
  return c_;
}

bool LatticeModel::CheckInput()
{
  return number_of_dimensions_ == 0 || number_of_directions_ == 0 ||
      number_of_rows_ == 0 || number_of_columns_ == 0;
}
