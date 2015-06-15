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
  , double dt
  , double diffusion_coefficient
  , double kinematic_viscosity
  , bool is_cd
  , bool is_ns)
  : number_of_dimensions_ {num_dimensions},
    number_of_discrete_velocities_ {num_discrete_velocities},
    number_of_rows_ {num_rows},
    number_of_columns_ {num_cols},
    space_step_ {dx},
    time_step_ {dt},
    diffusion_coefficient_ {diffusion_coefficient},
    kinematic_viscosity_ {kinematic_viscosity},
    is_cd_ {is_cd},
    is_ns_ {is_ns}

{
  if (input_parameter_check_value_ == 0) {
    throw std::runtime_error("Zero value in input parameters");
  }
  else {
    c_ = space_step_ / time_step_;
    cs_sqr_ = c_ * c_ / 3.0;
    std::size_t lattice_size = (num_rows) * (num_cols);
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
    // tau_cd_ formula from "A new scheme for source term in LBGK model for
    // convection–diffusion equation"
    if (is_cd) tau_cd_ = 0.5 + diffusion_coefficient / cs_sqr_ / dt;
    // tau_ns_ formula from "Discrete lattice effects on the forcing term in the
    // lattice Boltzmann method" Guo2002
    if (is_ns) tau_ns_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt;
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
  for (auto &lat : lattice) {
    if (initial_values.size() != lat.size()) {
      throw std::runtime_error("Depth mismatch");
    }
    else {
      lat = initial_values;
    }
  }  // lat
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
  if (src_position.size() != src_strength.size())
      throw std::runtime_error("Insufficient source information");
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
    lattice_src[src_pos[1] * nx + src_pos[0]] = *it_strength++;
  }  // src_pos
}

void Lattice::InitSrc(std::vector<std::vector<double>> &lattice_src
  , const std::vector<std::vector<unsigned>> &src_position
  , const std::vector<std::vector<double>> &src_strength)
{
  if (src_position.size() != src_strength.size())
      throw std::runtime_error("Insufficient source information");
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
    lattice_src[src_pos[1] * nx + src_pos[0]] = *it_strength++;
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
  for (auto n = 0u; n < nx * ny; ++n) {
    double u_sqr = u_[n][0] * u_[n][0] + u_[n][1] * u_[n][1];
    u_sqr /= 2. * cs_sqr;
    for (auto i = 0u; i < nc; ++i) {
      double c_dot_u = u_[n][0] * e_[i][0] + u_[n][1] * e_[i][1];
      c_dot_u /= cs_sqr / c_;
      lattice_eq[n][i] = omega_[i] * rho[n] * (1. + c_dot_u *
          (1. + c_dot_u / 2.) - u_sqr);
    }  // i
  }  // n
}

void Lattice::BoundaryCondition(const std::vector<std::vector<double>> &lattice
  , std::vector<std::vector<double>> &boundary)
{
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  std::vector<std::vector<double>> left_boundary(ny,
      std::vector<double>(9, 0.0));
  std::vector<std::vector<double>> right_boundary(ny,
      std::vector<double>(9, 0.0));
  std::vector<std::vector<double>> top_boundary(nx,
      std::vector<double>(9, 0.0));
  std::vector<std::vector<double>> bottom_boundary(nx,
      std::vector<double>(9, 0.0));
  std::vector<std::vector<double>> corner_boundary(4,
      std::vector<double>(9, 0.0));
  // Periodic boundary condition on left and right
  for (auto y = 0u; y < ny; ++y) {
    auto n = y * nx;
    left_boundary[y][E] = lattice[n + nx - 1][E];
    left_boundary[y][NE] = lattice[n + nx - 1][NE];
    left_boundary[y][SE] = lattice[n + nx - 1][SE];
    right_boundary[y][W] = lattice[n][W];
    right_boundary[y][NW] = lattice[n][NW];
    right_boundary[y][SW] = lattice[n][SW];
  }  // y
  // no-slip boundary condition on top and bottom
  for (auto x = 0u; x < nx; ++x) {
    auto n = (ny - 1) * (nx);
    top_boundary[x][S] = lattice[x + n][N];
    bottom_boundary[x][N] = lattice[x][S];
    if (x == 0) {
      top_boundary[x][SW] = lattice[nx * ny - 1][NE];
      bottom_boundary[x][NW] = lattice[nx - 1][SE];
    }
    else {
      top_boundary[x][SW] = lattice[x + n - 1][NE];
      bottom_boundary[x][NW] = lattice[x - 1][SE];
    }
    if (x == nx - 1) {
      top_boundary[x][SE] = lattice[n][NW];
      bottom_boundary[x][NE] = lattice[0][SW];
    }
    else {
      top_boundary[x][SE] = lattice[x + n + 1][NW];
      bottom_boundary[x][NE] = lattice[x + 1][SW];
    }
  }  // x
  // mixture of boundaries at the corners
  corner_boundary[0][NE] = lattice[0][SW];
  corner_boundary[1][NW] = lattice[nx - 1][SE];
  corner_boundary[2][SE] = lattice[(ny - 1) * nx][NW];
  corner_boundary[3][SW] = lattice[ny * nx - 1][NE];

  boundary.insert(end(boundary), begin(left_boundary), end(left_boundary));
  boundary.insert(end(boundary), begin(right_boundary), end(right_boundary));
  boundary.insert(end(boundary), begin(top_boundary), end(top_boundary));
  boundary.insert(end(boundary), begin(bottom_boundary), end(bottom_boundary));
  boundary.insert(end(boundary), begin(corner_boundary), end(corner_boundary));
}

std::vector<std::vector<double>> Lattice::Stream(
    const std::vector<std::vector<double>> &lattice
  , const std::vector<std::vector<double>> &boundary)
{
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  auto nc = GetNumberOfDiscreteVelocities();
  auto temp_lattice(lattice);
  std::vector<std::vector<double>> lattice_with_boundary((nx + 2) * (ny + 2),
      std::vector<double>(nc, 0.0));
  // assemble lattice and boundary into a temp lattice with boundary
  // surrounding it
  for (auto y = 1u; y < ny + 1; ++y) {
    // left and right boundary
    lattice_with_boundary[y * (nx + 2)] = boundary[y - 1];
    lattice_with_boundary[y * (nx + 2) + nx + 1] = boundary[ny + y - 1];
    // main lattice
    for (auto x = 1u; x < nx + 1; ++x) {
      auto n = y * (nx + 2) + x;
      auto m = (y - 1) * nx + x - 1;
      lattice_with_boundary[n] = lattice[m];
      if (y == 1) {
        // top and bottom boundary
        lattice_with_boundary[(ny + 1) * (nx + 2) + x] =
            boundary[2 * ny + x - 1];
        lattice_with_boundary[x] = boundary[2 * ny + nx + x - 1];
      }
    }  // x
  }  // y
  // corner boundary
  auto corner_index = 2 * ny + 2 * nx;
  lattice_with_boundary[0] = boundary[corner_index];
  lattice_with_boundary[nx + 1] = boundary[corner_index + 1];
  lattice_with_boundary[(ny + 1) * (nx + 2)] = boundary[corner_index + 2];
  lattice_with_boundary[(ny + 2) * (nx + 2) - 1] = boundary[corner_index + 3];
  for (auto y = 1u; y < ny + 1; ++y) {
    for (auto x = 1u; x < nx + 1; ++x) {
      auto n = y * (nx + 2) + x;
      auto m = (y - 1) * nx + x - 1;
      temp_lattice[m][E] = lattice_with_boundary[n - 1][E];
      temp_lattice[m][N] = lattice_with_boundary[n - (nx + 2)][N];
      temp_lattice[m][W] = lattice_with_boundary[n + 1][W];
      temp_lattice[m][S] = lattice_with_boundary[n + (nx + 2)][S];
      temp_lattice[m][NE] = lattice_with_boundary[n - (nx + 2) - 1][NE];
      temp_lattice[m][NW] = lattice_with_boundary[n - (nx + 2) + 1][NW];
      temp_lattice[m][SW] = lattice_with_boundary[n + (nx + 2) + 1][SW];
      temp_lattice[m][SE] = lattice_with_boundary[n + (nx + 2) - 1][SE];
    }  // x
  }  // y
  return temp_lattice;
}

void Lattice::Collide(std::vector<std::vector<double>> &lattice
  , const std::vector<std::vector<double>> &lattice_eq
  , std::vector<double> &src)
{
  if (lattice.size() != lattice_eq.size())
      throw std::runtime_error("Lattice size mismatch");
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  auto nc = GetNumberOfDiscreteVelocities();
//  auto src_i(lattice);
  for (auto n = 0u; n < nx * ny; ++n) {
    for (auto i = 0u; i < nc; ++i) {
      double c_dot_u = u_[n][0] * e_[i][0] + u_[n][1] * e_[n][1];
      c_dot_u /= cs_sqr_ / c_;
      // Source term using forward scheme, theta = 0
      auto src_i = omega_[i] * src[n] * (1.0 + (1.0 - 0.5 / tau_cd_) *
          c_dot_u);
      lattice[n][i] += (lattice_eq[n][i] - lattice[n][i]) / tau_cd_ +
          time_step_ * src_i;
    }  // i
    if (is_instant_) src[n] = 0.0;
  }  // n
}

void Lattice::Collide(std::vector<std::vector<double>> &lattice
  , const std::vector<std::vector<double>> &lattice_eq
  , std::vector<std::vector<double>> &src)
{
  if (lattice.size() != lattice_eq.size())
      throw std::runtime_error("Lattice size mismatch");
  auto cs_sqr = c_ * c_ / 3.0;

}

std::vector<double> Lattice::Flip(const std::vector<double> &lattice)
{
  auto flipped_lattice(lattice);
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  for (int y = ny - 1, y_flipped = 0; y > -1; --y, ++y_flipped) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      auto n_flipped = y_flipped * nx + x;
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
  for (int y = ny - 1, y_flipped = 0; y > -1; --y, ++y_flipped) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      auto n_flipped = y_flipped * nx + x;
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
    if (++counter % nx == 0) {
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
  // 0 for depth 9, 1 for depth 2, 2 for boundary, 3 for big
  switch (which_to_print) {
    case 0: {
      auto nc = GetNumberOfDiscreteVelocities();
      // row of lattice
      for (int y = ny - 1; y > -1; --y) {
        // rows in the Q9 square
        for (auto i = 0u; i < nc / 3; ++i) {
          // column of lattice
          for (auto x = 0u; x < nx; ++x) {
            int n = y * nx + x;
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
        if (++counter % nx == 0) {
          std::cout << std::endl;
          counter = 0;
        }
      }  // lat
      std::cout << std::endl;
      break;
    }
    case 2: {
      auto nc = GetNumberOfDiscreteVelocities();
      std::vector<unsigned> length = {ny, ny, nx, nx, 4};
      // row of lattice
      for (auto y = 0; y < 5; ++y) {
        // rows in the Q9 square
        for (auto i = 0u; i < nc / 3; ++i) {
          // column of lattice
          for (auto x = 0u; x < length[y]; ++x) {
            unsigned n;
            if (y == 0 || y == 1) {
              n = y * ny + x;
            }
            else if (y == 2 || y == 3) {
              n = 2 * ny + (y - 2) * nx + x;
            }
            else {
              n = 2 * ny + 2 * nx + x;
            }
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
    case 3: {
      auto nc = GetNumberOfDiscreteVelocities();
      // row of lattice
      for (int y = ny + 1; y > -1; --y) {
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
    default: {
      throw std::runtime_error("Not a 2D lattice");
      break;
    }
  }
}
