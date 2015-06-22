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
    auto nc = lm_.GetNumberOfDimensions();
    auto nd = lm_.GetNumberOfDirections();
    std::size_t lattice_size = ny * nx;
    std::vector<double> length_d(nd, 0.0);
    std::vector<double> length_q(nc, 0.0);
    // Initializing variables with initial values
    if (is_ns_) {
      f = ns_.lattice_eq;
      boundary_f.assign(2 * ny + 2 * nx + 4, length_q);
    }
    if (is_cd_) {
      g = cd_.lattice_eq;
      boundary_g.assign(2 * ny + 2 * nx + 4, length_q);
//      std::cout << lm.test_value << std::endl;
    }
    if (has_obstacles) {
      obstacles.assign(lattice_size, false);
      LatticeBoltzmann::Init(obstacles, obstacles_position);
    }
  }
}

std::vector<double> LatticeBoltzmann::GetRhoF() const
{
  return ns_.GetRho();
}

std::vector<double> LatticeBoltzmann::GetRhoG() const
{
  return cd_.GetRho();
}

std::vector<std::vector<double>> LatticeBoltzmann::GetVelocity() const
{
  return (is_ns_) ? ns_.GetVelocity() : cd_.GetVelocity();
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
