#include "LatticeD2Q9.hpp"
#include "LatticeModel.hpp"
LatticeD2Q9::LatticeD2Q9(std::size_t num_rows
  , std::size_t num_cols
  , double dx
  , double dt)
{
  // can't use initializer list because inheritance provides access but doesn't
  // create member variables in derived class, so have to reassign them values
  // the old fashioned way
  // https://stackoverflow.com/questions/6986798/subtle-c-inheritance-error-with
  // -protected-fields
  number_of_dimensions_ = 2;
  number_of_directions_ = 9;
  number_of_rows_ = num_rows;
  number_of_columns_ = num_cols;
  space_step_ = dx;
  time_step_ = dt;
  e = {{0, 0},
       {1, 0}, {0, 1}, {-1, 0}, {0, -1},
       {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  omega = {16.0 / 36.0,
           4.0 / 36.0, 4.0 / 36.0, 4.0 / 36.0, 4.0 / 36.0,
           1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
}
