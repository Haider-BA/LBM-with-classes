#ifndef WRITE_TO_CMGUI_HPP_
#define WRITE_TO_CMGUI_HPP_

#include <vector>

void WriteToCmgui(const std::vector<std::vector<double>> &f
  , const std::vector<std::vector<double>> &g
  , const std::vector<std::vector<double>> &u
  , std::size_t nx
  , std::size_t ny
  , int time
  , double rho0
  , double c
  , double cs_sqr);

#endif // WRITE_TO_CMGUI_HPP_
