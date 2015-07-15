#ifndef LATTICEBOLTZMANN_HPP_
#define LATTICEBOLTZMANN_HPP_
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeModel.hpp"

class LatticeBoltzmann {
 public:
  LatticeBoltzmann(double t_total
    , LatticeModel &lm
    , CollisionNS &ns);
  LatticeBoltzmann(double t_total
    , LatticeModel &lm
    , CollisionCD &cd);
  LatticeBoltzmann(double t_total
    , LatticeModel &lm
    , CollisionNS &ns
    , CollisionCD &cd);
  std::vector<std::vector<double>> Stream(
    const std::vector<std::vector<double>> &lattice);
 private:
  // 6  2  5  ^
  //  \ | /   |
  // 3--0--1  |
  //  / | \   |
  // 7  4  8  +------->
  enum Directions {
    E = 1,
    N,
    W,
    S,
    NE,
    NW,
    SW,
    SE
  };
  double total_time_;
  // LatticeModel to take care of dims, dirs, rows, cols and discrete e vectors
  // by reference, similar to by pointer
  // https://stackoverflow.com/questions/9285627/is-it-possible-to-pass-derived-
  // classes-by-reference-to-a-function-taking-base-cl
  LatticeModel &lm_;
  // Collision models to take care of eq, rho, source
  CollisionNS &ns_;
  CollisionCD &cd_;
  bool is_cd_;
  bool is_ns_;
};
#endif  // LATTICEBOLTZMANN_HPP_
