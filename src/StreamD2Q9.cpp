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
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto width = nx - 1;
  auto height = (ny - 1) * nx;
  auto temp_df = df;
  // Streaming
  for (auto n = 0u; n < nx * ny; ++n) {
    auto left = n % nx == 0;
    auto right = n % nx == nx - 1;
    auto bottom = n / nx == 0;
    auto top = n / nx == ny - 1;

    temp_df[n][E] = df[left ? n : n - 1][E];
    temp_df[n][N] = df[bottom ? n : n - nx][N];
    temp_df[n][W] = df[right ? n : n + 1][W];
    temp_df[n][S] = df[top ? n : n + nx][S];
    temp_df[n][NE] = df[(bottom || left)? n: n - nx - 1][NE];
    temp_df[n][NW] = df[(bottom || right) ? n : n -nx + 1][NW];
    temp_df[n][SW] = df[(top || right) ? n : n + nx + 1][SW];
    temp_df[n][SE] = df[(top || left) ? n : n + nx - 1][SE];
  }  // n
  return temp_df;
}
