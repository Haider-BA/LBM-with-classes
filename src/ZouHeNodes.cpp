#include "ZouHeNodes.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>
#include "ValueNode.hpp"

ZouHeNodes::ZouHeNodes(bool is_prestream
  , CollisionModel &cm
  , LatticeModel &lm)
  : BoundaryNodes(is_prestream, lm),
    nodes_ {}
{}

void ZouHeNodes::AddNode(std::size_t x
  , std::size_t y
  , double u_lid
  , double v_lid)
{
  // should just try to implement normally using switch statement first
  auto nx = lm_.GetNumberOfColumns();
  auto ny = lm_.GetNumberOfRows();
  auto left = x == 0;
  auto right = x == nx - 1;
  auto bottom = y == 0;
  auto top = y == ny - 1;
//  knowns_ = {!((left || right) && (top || bottom)),
//      right, top, left, bottom,
//      right || top, left || top, left || bottom, right || bottom};
  auto n = y * nx + x;
  auto side = -1;
  if (right) side = 0;
  if (top) side = 1;
  if (left) side = 2;
  if (bottom) side = 3;
  nodes_.push_back(ValueNode(x, y, nx, u_lid, v_lid, ((top || bottom) &&
      (left || right)) ? true : false, side));
}

void ZouHeNodes::UpdateNodes(std::vector<std::vector<double>> &df)
{
  for (auto n = 0u; n < nodes_.size(); ++n) {
    if (nodes_[n].b1) {
      ZouHeNodes::UpdateCorner(df, nodes_[n].n);
    }
    else {
      ZouHeNodes::UpdateSide(df, nodes_[n]);
    }
  }
}

void ZouHeNodes::UpdateSide(std::vector<std::vector<double>> &df
  , ValueNode &node)
{
  auto n = node.n;
  auto vel = node.v1;
  switch(node.i1) {
    case 0: {
      auto node_rho = -1.0 * vel[0] + df[n][0] + df[n][N] + df[n][S] + 2.0 *
          (df[n][E] + df[n][NE] + df[n][SE]);
      auto df_diff = 0.5 * (df[n][S] - df[n][N]);
      df[n][W] = df[n][E] - 2.0 / 3.0 * vel[0];
      df[n][NW] = df[n][SE] + df_diff - vel[0] / 6.0 + vel[1] / 2.0;
      df[n][SW] = df[n][NE] - df_diff - vel[0] / 6.0 - vel[1] / 2.0;
      break;
    }
    case 1: {
      auto node_rho = -1.0 * vel[1] + df[n][0] + df[n][E] + df[n][W] + 2.0 *
          (df[n][N] + df[n][NE] + df[n][NW]);
      auto df_diff = 0.5 * (df[n][E] - df[n][W]);
      df[n][S] = df[n][N] - 2.0 / 3.0 * vel[1];
      df[n][SW] = df[n][NE] + df_diff - vel[0] / 2.0 - vel[1] / 6.0;
      df[n][SE] = df[n][NW] - df_diff + vel[0] / 2.0 - vel[1] / 6.0;
      break;
    }
    case 2: {
      auto node_rho = vel[0] + df[n][0] + df[n][N] + df[n][S] + 2.0 *
          (df[n][W] + df[n][NW] + df[n][SW]);
      auto df_diff = 0.5 * (df[n][S] - df[n][N]);
      df[n][E] = df[n][W] + 2.0 / 3.0 * vel[0];
      df[n][NE] = df[n][SW] + df_diff + vel[0] / 6.0 + vel[1] / 2.0;
      df[n][SE] = df[n][NW] - df_diff + vel[0] / 6.0 - vel[1] / 2.0;
      break;
    }
    case 3: {
      auto node_rho = vel[1] + df[n][0] + df[n][E] + df[n][W] + 2.0 *
          (df[n][S] + df[n][SW] + df[n][SE]);
      auto df_diff = 0.5 * (df[n][W] - df[n][E]);
      df[n][N] = df[n][S] + 2.0 / 3.0 * vel[1];
      df[n][NE] = df[n][SW] + df_diff + vel[0] / 6.0 + vel[1] / 2.0;
      df[n][NW] = df[n][SE] - df_diff - vel[0] / 6.0 + vel[1] / 2.0;
      break;
    }
    default: {
      throw std::runtime_error("Not a side");
    }
  }
//  auto n = node.n;
//  // v1 is vector for lid velocity
//  // i1 is int for which side, 1 = right, 2 = top, 3 = left, 4 = bottom
//  auto rho = node.v1[node.i1 % 2] + df[n][0] +
//      (node.i1 % 2 ? df[n][1] + df[n][3] : df[n][2] + df[n][4]) +
//      2.0 * (df[n][node.i1] + df[n][node.i1 + (node.i1 == 1) ? 7 : 3] +
//      df[n][node.i1 + 4]);
//  df[n][node.i1] =

}

void ZouHeNodes::UpdateCorner(std::vector<std::vector<double>> &df
    , std::size_t n)
{
  throw std::runtime_error("corner not implemented");
}
