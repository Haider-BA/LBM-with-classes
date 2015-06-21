#include "Collision.hpp"
#include <iostream>
#include "LatticeModel.hpp"

Collision::Collision(LatticeModel &lm)
  : lm_ (lm)
{
  auto dx = lm.GetSpaceStep();
  auto dt = lm.GetTimeStep();
  c_ = dx / dt;
  cs_sqr_ = c_ * c_ / 3.0;
}

void Collision::Collide()
{
  std::cout << lattice_eq_[0][0] << std::endl;
}
