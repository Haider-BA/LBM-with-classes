#ifndef PRINTING_HPP_
#define PRINTING_HPP_
#include <iomanip>  // std::setprecision
#include <iostream>
#include <stdexcept>  // std::runtime_error
#include <vector>

/**
 * Flips the lattice into the proper orientation for printing
 * \param lattice lattice stored row-wise
 * \param nx number of columns of the lattice
 * \return flipped lattice
 */
template <typename T, typename U>
T Flip(const T &lattice
  , U nx)
{
  auto flipped_lattice(lattice);
  auto ny = lattice.size() / nx;
  for (auto n = 0u; n < nx * ny; ++n) {
    auto n_flipped = (ny - 1 - n / nx) * nx + n % nx;
    flipped_lattice[n_flipped] = lattice[n];
  }  // n
  return flipped_lattice;
}

/**
 * Reorganize lattice using new index format for easier printing
 * 0  1  2
 *  \ | /
 * 3--4--5
 *  / | \
 * 6  7  8
 * \param lattice lattice stored row-wise
 * \return reorganized lattice
 */
template <typename T>
T Organize(const T &lattice)
{
  auto organized_lattice(lattice);
  auto index = 0;
  for (auto &org : organized_lattice) {
    // 3 and 8 are in organized position
    org[0] = lattice[index][6];
    org[1] = lattice[index][2];
    org[2] = lattice[index][5];
    org[4] = lattice[index][0];
    org[5] = lattice[index][1];
    org[6] = lattice[index][7];
    org[7] = lattice[index][4];
    ++index;
  }  // org
  return organized_lattice;
}

/**
 * Prints a depth 1 lattice
 * \param lattice depth 1 lattice stored row-wise
 * \param nx number of columns of the lattice
 */
template <typename T, typename U>
void Print(const T &lattice
  , U nx)
{
  auto counter = 0;
  auto flipped_lattice = Flip(lattice, nx);
  for (auto node : flipped_lattice) {
    std::cout << std::fixed << std::setprecision(2) << node << " ";
    if (++counter % nx == 0) {
      std::cout << std::endl;
      counter = 0;
    }
  }  // lat
  std::cout << std::endl;
}

/**
 * Print a multi-depth lattice
 * \param lattice stored row-wise
 * \param nx number of columns of lattice
 * \param ny number of rows of lattice
 */
template <typename T, typename U>
void Print(const T &lattice
  , U nx
  , U ny)
{
  auto depth = lattice.at(0).size();
  auto code = depth + lattice.size() % (nx * ny);
  // 2 for depth 2, 9 for depth 9, default for boundary
  switch (code) {
    case 2: {
      int counter = 0;
      auto flipped_lattice = Flip(lattice, nx);
      for (auto node : flipped_lattice) {
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
    case 9: {
      auto lat = Flip(Organize(lattice), nx);
      // row of lattice
      for (auto y = 0u; y < ny; ++y) {
        // rows in the Q9 square
        for (auto i = 0u; i < depth / 3; ++i) {
          // column of lattice
          for (auto x = 0u; x < nx; ++x) {
            auto n = y * nx + x;
            auto index = i * 3;
            std::cout << std::fixed << std::setprecision(2)
                      << lat[n][index] << " "
                      << lat[n][index + 1] << " "
                      << lat[n][index + 2] << "  ";
          }  // x
          std::cout << std::endl;
        }  // i
        std::cout << std::endl;
      }  // y
      std::cout << std::endl;
      break;
    }
    default: {
      std::vector<U> length = {ny, ny, nx, nx, 4};
      auto lat = Organize(lattice);
      // row of boundary, print 4 lines, left, right, top, bottom and corners
      for (auto y = 0u; y < 5; ++y) {
        // rows in the Q9 square
        for (auto i = 0u; i < depth / 3; ++i) {
          // column of lattice
          for (auto x = 0u; x < length[y]; ++x) {
            auto n = 2 * ny + 2 * nx + x;
            auto index = i * 3;
            if (y == 0 || y == 1) n = y * ny + x;
            if (y == 2 || y == 3) n = 2 * ny + (y - 2) * nx + x;
            std::cout << std::fixed << std::setprecision(2)
                      << lat[n][index] << " "
                      << lat[n][index + 1] << " "
                      << lat[n][index + 2] << "  ";
          }  // x
          std::cout << std::endl;
        }  // i
        std::cout << std::endl;
      }  // y
      std::cout << std::endl;
      break;
    }
  }
}
#endif  // PRINTING_HPP_
