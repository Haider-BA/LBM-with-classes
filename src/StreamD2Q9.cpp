#include "StreamD2Q9.hpp"
#include <iostream>
#include <vector>
#include "LatticeModel.hpp"
#include "StreamModel.hpp"

StreamD2Q9::StreamD2Q9(LatticeModel &lm)
  : StreamModel(lm)
{}

std::vector<std::vector<double>> StreamD2Q9::Stream(
    const std::vector<std::vector<double>> &df)
{
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  auto temp_df = df;
  // Streaming
  for (auto n = 0u; n < nx * ny; ++n) {
    const auto left = n % nx == 0;
    const auto right = n % nx == nx - 1;
    const auto bottom = n / nx == 0;
    const auto top = n / nx == ny - 1;
    if (!left) temp_df[n][E] = df[n - 1][E];
    if (!bottom) temp_df[n][N] = df[n - nx][N];
    if (!right) temp_df[n][W] = df[n + 1][W];
    if (!top) temp_df[n][S] = df[n + nx][S];
    if (!(bottom || left)) temp_df[n][NE] = df[n - nx - 1][NE];
    if (!(bottom || right)) temp_df[n][NW] = df[n - nx + 1][NW];
    if (!(top || right)) temp_df[n][SW] = df[n + nx + 1][SW];
    if (!(top || left)) temp_df[n][SE] = df[n + nx - 1][SE];
  }  // n
  return temp_df;
}
