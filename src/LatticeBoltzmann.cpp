#include "LatticeBoltzmann.hpp"
#include <cmath>  // std::fmod
#include <iomanip>  // std::setprecision
#include <iostream>
#include <stdexcept>  // std::runtime_error
#include <vector>
#include "Algorithm.hpp"
#include "CollisionModel.hpp"
#include "LatticeModel.hpp"
#include "Printing.hpp"
#include "WriteResultsCmgui.hpp"

// cant use braces to initialize reference cuz gcc bug
// https://stackoverflow.com/questions/10509603/why-cant-i-initialize-a-
// reference-in-an-initializer-list-with-uniform-initializ
LatticeBoltzmann::LatticeBoltzmann(LatticeModel &lm
  , CollisionModel &cm
  , StreamModel &sm)
  : df {cm.edf},
    lm_ (lm),
    cm_ (cm),
    sm_ (sm),
    bn_ {}
{}

void LatticeBoltzmann::AddBoundaryNodes(BoundaryNodes *bn)
{
  bn_.push_back(bn);
}

void LatticeBoltzmann::TakeStep()
{
  cm_.ComputeMacroscopicProperties(df);
  cm_.ComputeEq();
  cm_.Collide(df);
  for (auto bdr : bn_) {
    if (bdr->prestream) bdr->UpdateNodes(df, false);
  }  // bdr
  df = sm_.Stream(df);
  for (auto bdr : bn_) {
    if (bdr-> during_stream) bdr->UpdateNodes(df, true);
    if (!bdr->prestream) bdr->UpdateNodes(df, false);
  }  // bdr
}
