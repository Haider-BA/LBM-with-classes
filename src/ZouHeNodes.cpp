#include "ZouHeNodes.hpp"
#include <iostream>
#include <vector>
#include "ValueNode.hpp"

ZouHeNodes::ZouHeNodes(bool is_prestream
  , CollisionModel &cm
  , LatticeModel &lm)
  : BoundaryNodes(is_prestream, lm),
    is_corner_ {}
{
  ;
}

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
  knowns_ = {!((left || right) && (top || bottom)),
      right, top, left, bottom,
      right || top, left || top, left || bottom, right || bottom};
  auto n = y * nx + x;
  if (left) {
    nodes_.push_back(ValueNode(x, y, nx, u_lid, v_lid, ((top || bottom)) ?
        true : false));
  }
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
  ;
}
