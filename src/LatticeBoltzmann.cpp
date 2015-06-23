#include "LatticeBoltzmann.hpp"
#include <iomanip> // std::setprecision
#include <iostream>
#include <stdexcept>  // std::runtime_error
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeModel.hpp"

//LatticeBoltzmann::LatticeBoltzmann()
//{
//  throw std::runtime_error("Declaring uninitialized LBM is not allowed");
//}

// cant use braces to initialize reference cuz gcc bug
// https://stackoverflow.com/questions/10509603/why-cant-i-initialize-a-
// reference-in-an-initializer-list-with-uniform-initializ
LatticeBoltzmann::LatticeBoltzmann(double t_total
  , const std::vector<std::vector<std::size_t>> &obstacles_position
  , bool is_ns
  , bool is_cd
  , bool is_instant
  , bool has_obstacles
  , LatticeModel &lm
  , CollisionNS &ns
  , CollisionCD &cd)
  : total_time_ {t_total},
    is_ns_ {is_ns},
    is_cd_ {is_cd},
    is_instant_ {is_instant},
    has_obstacles_ {has_obstacles},
    lm_ (lm),
    ns_ (ns),
    cd_ (cd)
{
  if (CheckParameters()) {
    throw std::runtime_error("Zero value in input parameters");
  }
  else {
    // Initializing all variables with zero values
    auto nx = lm_.GetNumberOfColumns();
    auto ny = lm_.GetNumberOfRows();
    auto nc = lm_.GetNumberOfDirections();
    auto nd = lm_.GetNumberOfDimensions();
    std::size_t lattice_size = ny * nx;
    std::vector<double> length_d(nd, 0.0);
    std::vector<double> length_q(nc, 0.0);
    // Initializing variables with initial values
    if (is_ns_) {
      f = ns_.lattice_eq;
//      boundary_f.assign(2 * ny + 2 * nx + 4, length_q);
    }
    if (is_cd_) {
      g = cd_.lattice_eq;
//      boundary_g.assign(2 * ny + 2 * nx + 4, length_q);
//      std::cout << lm.test_value << std::endl;
    }
    if (has_obstacles) {
      obstacles.assign(lattice_size, false);
      LatticeBoltzmann::Init(obstacles, obstacles_position);
    }
  }
}

void LatticeBoltzmann::Init(std::vector<bool> &lattice
  , const std::vector<std::vector<std::size_t>> &position)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nd = lm_.GetNumberOfDimensions();
  if(position[0].size() != nd)
      throw std::runtime_error("Insufficient position information");
  for (auto pos : position) {
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    lattice[pos[1] * nx + pos[0]] = true;
  }  // pos
}

std::vector<std::vector<double>> LatticeBoltzmann::BoundaryCondition(
    const std::vector<std::vector<double>> &lattice)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nc = lm_.GetNumberOfDirections();
  std::vector<double> length_q(nc, 0.0);
  std::vector<std::vector<double>> left_boundary(ny, length_q);
  std::vector<std::vector<double>> top_boundary(nx, length_q);
  std::vector<std::vector<double>> corner_boundary(4, length_q);
  // initialize some boundaries to mirror their counterparts
  auto right_boundary(left_boundary);
  auto bottom_boundary(top_boundary);
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
    auto n = (ny - 1) * nx;
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

  // append the boundaries of different edges to the main boundary vector
  auto boundary(left_boundary);
  boundary.insert(end(boundary), begin(right_boundary), end(right_boundary));
  boundary.insert(end(boundary), begin(top_boundary), end(top_boundary));
  boundary.insert(end(boundary), begin(bottom_boundary), end(bottom_boundary));
  boundary.insert(end(boundary), begin(corner_boundary), end(corner_boundary));
  return boundary;
}

std::vector<std::vector<double>> LatticeBoltzmann::Stream(
    const std::vector<std::vector<double>> &lattice
  , const std::vector<std::vector<double>> &boundary)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto nc = lm_.GetNumberOfDirections();
  auto temp_lattice(lattice);
  std::vector<std::vector<double>> lattice_with_boundary((nx + 2) * (ny + 2),
      std::vector<double>(nc, 0.0));
  // assemble lattice and boundary into a lattice with boundary
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

bool LatticeBoltzmann::CheckParameters()
{
  return total_time_ < 1e-20 || !(is_ns_ || is_cd_);
}

std::vector<double> LatticeBoltzmann::Flip(const std::vector<double> &lattice)
{
  auto flipped_lattice(lattice);
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (int y = ny - 1, y_flipped = 0; y > -1; --y, ++y_flipped) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      auto n_flipped = y_flipped * nx + x;
      flipped_lattice[n_flipped] = lattice[n];
    }  // x
  }  // y
  return flipped_lattice;
}

std::vector<bool> LatticeBoltzmann::Flip(const std::vector<bool> &lattice)
{
  auto flipped_lattice(lattice);
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (int y = ny - 1, y_flipped = 0; y > -1; --y, ++y_flipped) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      auto n_flipped = y_flipped * nx + x;
      flipped_lattice[n_flipped] = lattice[n];
    }  // x
  }  // y
  return flipped_lattice;
}

std::vector<std::vector<double>> LatticeBoltzmann::Flip(
    const std::vector<std::vector<double>> &lattice)
{
  auto flipped_lattice(lattice);
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  for (int y = ny - 1, y_flipped = 0; y > -1; --y, ++y_flipped) {
    for (auto x = 0u; x < nx; ++x) {
      auto n = y * nx + x;
      auto n_flipped = y_flipped * nx + x;
      flipped_lattice[n_flipped] = lattice[n];
    }  // x
  }  // y
  return flipped_lattice;
}

void LatticeBoltzmann::Print(const std::vector<double> &lattice)
{
  auto nx = lm_.GetNumberOfColumns();
  int counter = 0;
  auto flipped_lattice = LatticeBoltzmann::Flip(lattice);
  for (auto node : flipped_lattice) {
    std::cout << std::fixed << std::setprecision(2) << node << " ";
    if (++counter % nx == 0) {
      std::cout << std::endl;
      counter = 0;
    }
  }  // lat
  std::cout << std::endl;
}

void LatticeBoltzmann::Print(const std::vector<bool> &lattice)
{
  auto nx = lm_.GetNumberOfColumns();
  int counter = 0;
  auto flipped_lattice = LatticeBoltzmann::Flip(lattice);
  for (auto node : flipped_lattice) {
    std::cout << std::fixed << std::setprecision(2) << node << " ";
    if (++counter % nx == 0) {
      std::cout << std::endl;
      counter = 0;
    }
  }  // lat
  std::cout << std::endl;
}

void LatticeBoltzmann::Print(int which_to_print
  , const std::vector<std::vector<double>> &lattice)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  // 0 for depth 9, 1 for depth 2, 2 for boundary, 3 for big
  switch (which_to_print) {
    case 0: {
      auto nc = lm_.GetNumberOfDirections();
      // row of lattice
      for (int y = ny - 1; y > -1; --y) {
        // rows in the Q9 square
        for (auto i = 0u; i < nc / 3; ++i) {
          // column of lattice
          for (auto x = 0u; x < nx; ++x) {
            int n = y * nx + x;
            if (lattice[n].size() != nc)
                throw std::runtime_error("Wrong depth");
            if (i == 0) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][6] << " "
                      << lattice[n][2] << " "
                      << lattice[n][5] << " ";
            }
            else if (i == 1) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][3] << " "
                      << lattice[n][0] << " "
                      << lattice[n][1] << " ";
            }
            else if (i == 2) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][7] << " "
                      << lattice[n][4] << " "
                      << lattice[n][8] << " ";
            }
            std::cout << "  ";
          }  // x
          std::cout << std::endl;
        }  // i
        std::cout << std::endl;
      }  // y
      std::cout << std::endl;
      break;
    }
    case 1: {
      auto nd = lm_.GetNumberOfDimensions();
      int counter = 0;
      auto flipped_lattice = LatticeBoltzmann::Flip(lattice);
      for (auto node : flipped_lattice) {
        if (node.size() != nd) throw std::runtime_error("Wrong depth");
        std::cout << std::fixed << std::setprecision(2)
                  << node[0] << " " << node[1] << "  ";
        if (++counter % nx == 0) {
          std::cout << std::endl;
          counter = 0;
        }
      }  // lat
      std::cout << std::endl;
      break;
    }
    case 2: {
      auto nc = lm_.GetNumberOfDirections();
      std::vector<std::size_t> length = {ny, ny, nx, nx, 4};
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
            if (i == 0) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][6] << " "
                      << lattice[n][2] << " "
                      << lattice[n][5] << " ";
            }
            else if (i == 1) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][3] << " "
                      << lattice[n][0] << " "
                      << lattice[n][1] << " ";
            }
            else if (i == 2) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][7] << " "
                      << lattice[n][4] << " "
                      << lattice[n][8] << " ";
            }
            std::cout << "  ";
          }  // x
          std::cout << std::endl;
        }  // i
        std::cout << std::endl;
      }  // y
      std::cout << std::endl;
      break;
    }
    case 3: {
      auto nc = lm_.GetNumberOfDirections();
      // row of lattice
      for (int y = ny + 1; y > -1; --y) {
        // rows in the Q9 square
        for (auto i = 0u; i < nc / 3; ++i) {
          // column of lattice
          for (auto x = 0u; x < nx + 2; ++x) {
            int n = y * (nx + 2) + x;
            if (lattice[n].size() != nc)
                throw std::runtime_error("Wrong depth");
            if (i == 0) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][6] << " "
                      << lattice[n][2] << " "
                      << lattice[n][5] << " ";
            }
            else if (i == 1) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][3] << " "
                      << lattice[n][0] << " "
                      << lattice[n][1] << " ";
            }
            else if (i == 2) {
            std::cout << std::fixed << std::setprecision(2)
                      << lattice[n][7] << " "
                      << lattice[n][4] << " "
                      << lattice[n][8] << " ";
            }
            std::cout << "  ";
          }  // x
          std::cout << std::endl;
        }  // i
        std::cout << std::endl;
      }  // y
      std::cout << std::endl;
      break;
    }
    default: {
      throw std::runtime_error("Not a 2D lattice");
      break;
    }
  }
}
