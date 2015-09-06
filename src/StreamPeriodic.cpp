#include "StreamPeriodic.hpp"
#include <iostream>
#include <vector>
#include "LatticeModel.hpp"
#include "StreamModel.hpp"

StreamPeriodic::StreamPeriodic(LatticeModel &lm)
  : StreamModel(lm)
{}

std::vector<std::vector<double>> StreamPeriodic::Stream(
    const std::vector<std::vector<double>> &df)
{
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  const auto width = nx - 1;
  const auto height = (ny - 1) * nx;
  auto temp_df = df;
  // Streaming
  for (auto n = 0u; n < nx * ny; ++n) {
    auto left = n % nx == 0;
    auto right = n % nx == nx - 1;
    auto bottom = n / nx == 0;
    auto top = n / nx == ny - 1;
    temp_df[n][E] = df[left ? n + width : n - 1][E];
    temp_df[n][N] = df[bottom ? n + height : n - nx][N];
    temp_df[n][W] = df[right ? n - width : n + 1][W];
    temp_df[n][S] = df[top ? n - height : n + nx][S];
    if (left) {
      temp_df[n][NE] = df[bottom ? n + width + height : n + width - nx][NE];
      temp_df[n][SE] = df[top ? n + width - height : n + width + nx][SE];
    }
    else {
      temp_df[n][NE] = df[bottom ? n - 1 + height : n - 1 - nx][NE];
      temp_df[n][SE] = df[top ? n - 1 - height : n - 1 + nx][SE];
    }
    if (right) {
      temp_df[n][NW] = df[bottom ? n - width + height : n - width - nx][NW];
      temp_df[n][SW] = df[top ? n - width - height : n - width + nx][SW];
    }
    else {
      temp_df[n][NW] = df[bottom ? n + 1 + height : n + 1 - nx][NW];
      temp_df[n][SW] = df[top ? n + 1 - height : n + 1 + nx][SW];
    }
  }  // n
  return temp_df;
}
