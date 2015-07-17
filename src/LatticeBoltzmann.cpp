#include "LatticeBoltzmann.hpp"
#include <cstddef>  // NULL
#include <iostream>
#include <vector>
#include "CollisionModel.hpp"
#include "LatticeModel.hpp"
#include "LatticeD2Q9.hpp"
#include "Printing.hpp"

LatticeBoltzmann::LatticeBoltzmann(double t_total
  , LatticeModel &lm
  , CollisionModel &cm)
  : f {},
    g {},
    total_time_ {t_total},
    lm_ (lm),
    cm_ (cm),
    is_ns_ {cm.is_ns},
    is_cd_ {cm.is_cd}
{
  if (is_ns_) {
    auto nx = lm_.GetNumberOfColumns();
    auto ny = lm_.GetNumberOfRows();
    std::cout << "yes ns" << std::endl;
    f = cm.f_eq;
//    std::cout << cm.f_eq[1][5];
//    Print(cm_.f_eq, nx, ny);
  }
}

std::vector<std::vector<double>> LatticeBoltzmann::Stream(
    const std::vector<std::vector<double>> &lattice)
{
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto temp_lattice = lattice;
  for (auto n = 0u; n < nx * ny; ++n) {
    auto left = n % nx == 0;
    auto right = n % nx == nx - 1;
    auto top = n / nx == ny - 1;
    auto bottom = n / nx == 0;
    temp_lattice[n][E] = lattice[left ? n : n - 1][E];
    temp_lattice[n][N] = lattice[bottom ? n : n - nx][N];
    temp_lattice[n][W] = lattice[right ? n : n + 1][W];
    temp_lattice[n][S] = lattice[top ? n : n + nx][S];
    temp_lattice[n][NE] = lattice[(left || bottom) ? n : n - nx - 1][NE];
    temp_lattice[n][NW] = lattice[(right || bottom) ? n : n - nx + 1][NW];
    temp_lattice[n][SW] = lattice[(right || top) ? n : n + nx + 1][SW];
    temp_lattice[n][SE] = lattice[(left || top) ? n : n + nx - 1][SE];
  }  // n
  return temp_lattice;
}
