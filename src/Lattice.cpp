#include "Lattice.hpp"
#include <iomanip> // std::setprecision
#include <iostream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

Lattice::Lattice()
  : number_of_dimensions_ {0},
    number_of_discrete_velocities_ {0},
    number_of_rows_ {0},
    number_of_columns_ {0}
{}

Lattice::Lattice(std::size_t num_dimensions
  , std::size_t num_discrete_velocities
  , std::size_t num_rows
  , std::size_t num_cols
  , double dx
  , double dt)
  : number_of_dimensions_ {num_dimensions},
    number_of_discrete_velocities_ {num_discrete_velocities},
    number_of_rows_ {num_rows},
    number_of_columns_ {num_cols},
    space_step_ {dx},
    time_step_ {dt}

{
  if (input_parameter_check_value_ == 0) {
    throw std::runtime_error("Zero value in input parameters");
  }
  else {
    c_ = space_step_ / time_step_;
    std::size_t lattice_size = (num_rows + 2) * (num_cols + 2);
    std::vector<double> length_d(num_dimensions, 0.0);
    std::vector<double> length_q(num_discrete_velocities, 0.0);
    f_.assign(lattice_size, length_q);
    g_.assign(lattice_size, length_q);
    f_eq_.assign(lattice_size, length_q);
    g_eq_.assign(lattice_size, length_q);
    src_f_.assign(lattice_size, length_d);
    src_g_.assign(lattice_size, 0.0);
    rho_f_.assign(lattice_size, 0.0);
    rho_g_.assign(lattice_size, 0.0);
    u_.assign(lattice_size, length_d);
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

void Lattice::Init(std::vector<double> &lattice
    , double initial_value)
{
  for (auto &lat : lattice) lat = initial_value;
}

void Lattice::Init(std::vector<std::vector<double>> &lattice
    , const std::vector<double> &initial_values)
{
  if (initial_values.size() != number_of_dimensions_)
      throw std::runtime_error("Depth mismatch");
  for (auto &lat : lattice) lat = initial_values;
}

void Lattice::Init(std::vector<std::vector<double>> &lattice
  , const std::vector<std::vector<double>> &initial_lattice)
{
  if (initial_lattice.size() != lattice.size())
      throw std::runtime_error("Size mismatch");
  auto nc = GetNumberOfDiscreteVelocities();
  for (auto init_lat : initial_lattice) {
    if (init_lat.size() != nc) throw std::runtime_error("Depth mismatch");
  }  // init_lat
  lattice = initial_lattice;
}

void Lattice::InitSrc(std::vector<double> &lattice_src
    , const std::vector<std::vector<unsigned>> &src_position
    , const std::vector<double> &src_strength)
{
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  auto nd = GetNumberOfDimensions();
  auto it_strength = begin(src_strength);
  for (auto src_pos : src_position) {
    if(src_pos.size() != nd) {
      throw std::runtime_error("Insufficient position information");
    }
    else if (src_pos[0] > nx - 1) {
      throw std::runtime_error("x value out of range");
    }
    else if (src_pos[1] > ny - 1) {
      throw std::runtime_error("y value out of range");
    }
    auto n = src_pos[1] * (nx + 2) + src_pos[0] + nx + 3;
    lattice_src[n] = *it_strength++;
  }  // src_pos
}
void Lattice::ComputeEq(std::vector<std::vector<double>> &lattice_eq
  , const std::vector<double> &rho)
{
  auto nc = GetNumberOfDiscreteVelocities();
  for (auto lat : lattice_eq) {
    if (lat.size() != nc) throw std::runtime_error("Not depth 9 lattice.");
  }  // lat
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  double cs_sqr = c_ * c_ / 3.;
  for (auto y = 1u; y < ny + 1; ++y) {
    for (auto x = 1u; x < nx + 1; ++x) {
      auto n = y * (nx + 2) + x;
      double u_sqr = u_[n][0] * u_[n][0] + u_[n][1] * u_[n][1];
      u_sqr /= 2. * cs_sqr;
      for (auto i = 0u; i < nc; ++i) {
        double c_dot_u = u_[n][0] * e_[i][0] + u_[n][1] * e_[i][1];
        c_dot_u /= cs_sqr / c_;
        lattice_eq[n][i] = omega_[i] * rho[n] * (1. + c_dot_u *
            (1. + c_dot_u / 2.) - u_sqr);
      }  // i
    }  // x
  }  // y
}

std::vector<double> Lattice::Flip(const std::vector<double> &lattice)
{
  auto flipped_lattice(lattice);
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  for (int y = ny + 1, y_flipped = 0; y > -1; --y, ++y_flipped) {
    for (auto x = 0u; x < nx + 2; ++x) {
      auto n = y * (nx + 2) + x;
      auto n_flipped = y_flipped * (nx + 2) + x;
      flipped_lattice[n_flipped] = lattice[n];
    }
  }
  return flipped_lattice;
}

std::vector<std::vector<double>> Lattice::Flip(
    const std::vector<std::vector<double>> &lattice)
{
  auto flipped_lattice(lattice);
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  for (int y = ny + 1, y_flipped = 0; y > -1; --y, ++y_flipped) {
    for (auto x = 0u; x < nx + 2; ++x) {
      auto n = y * (nx + 2) + x;
      auto n_flipped = y_flipped * (nx + 2) + x;
      flipped_lattice[n_flipped] = lattice[n];
    }
  }
  return flipped_lattice;
}

void Lattice::Print(const std::vector<double> &lattice)
{
  auto nx = GetNumberOfColumns();
  int counter = 0;
  auto flipped_lattice = Lattice::Flip(lattice);
  for (auto lat : flipped_lattice) {
    std::cout << std::fixed << std::setprecision(2) << lat << " ";
    if (++counter % (nx + 2) == 0) {
      std::cout << std::endl;
      counter = 0;
    }
  }  // lat
  std::cout << std::endl;
}

void Lattice::Print(int which_to_print
  , const std::vector<std::vector<double>> &lattice)
{
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  // 0 for depth 9, 1 for depth 2
  switch (which_to_print) {
    case 0: {
      auto nc = GetNumberOfDiscreteVelocities();
      // row of lattice
      for (int y = ny + 1; y >-1; --y) {
        // rows in the Q9 square
        for (auto i = 0u; i < nc / 3; ++i) {
          // column of lattice
          for (auto x = 0u; x < nx + 2; ++x) {
            int n = y * (nx + 2) + x;
            if (lattice[n].size() != nc)
                throw std::runtime_error("Wrong depth");
            if (i == 0){
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][6] << " "
                      << lattice[n][2] << " "
                      << lattice[n][5] << " ";
            }
            else if (i == 1){
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][3] << " "
                      << lattice[n][0] << " "
                      << lattice[n][1] << " ";
            }
            else if (i == 2){
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][7] << " "
                      << lattice[n][4] << " "
                      << lattice[n][8] << " ";
            }
            std::cout << "  ";
          } // x
          std::cout << std::endl;
        } // i
        std::cout << std::endl;
      } // y
      std::cout << std::endl;
      break;
    }
    case 1: {
      auto nd = GetNumberOfDimensions();
      int counter = 0;
      auto flipped_lattice = Lattice::Flip(lattice);
      for (auto lat : flipped_lattice) {
        if (lat.size() != nd) throw std::runtime_error("Wrong depth");
        std::cout << std::fixed << std::setprecision(2)
                  << lat[0] << " " << lat[1] << "  ";
        if (++counter % (nx + 2) == 0) {
          std::cout << std::endl;
          counter = 0;
        }
      }  // lat
      std::cout << std::endl;
      break;
    }
    default: {
      throw std::runtime_error("Not a 2D lattice");
      break;
    }
  }
}
