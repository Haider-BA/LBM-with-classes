#include "LatticeBoltzmann.hpp"
#include <cstddef>  // NULL
#include <iostream>
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeModel.hpp"
#include "LatticeD2Q9.hpp"

LatticeBoltzmann::LatticeBoltzmann(double t_total
  , LatticeModel &lm
  , CollisionNS &ns
  , CollisionCD &cd)
  : total_time_ {t_total},
    lm_ (lm),
    ns_ (ns),
    cd_ (cd),
    is_cd_ {cd.is_implemented},
    is_ns_ {ns.is_implemented}
{
  if (is_cd_) std::cout << t_total << std::endl;
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
