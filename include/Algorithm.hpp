#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_
#include <cmath>  // fabs
#include <iostream>
#include <stdexcept>
#include <typeinfo>
#include <vector>

/**
 * Does dot product between 2 vectors of equal length
 * \param a first vector
 * \param b second vector
 * \return dot product of a and b
 */
template <typename T>
T InnerProduct(const std::vector<T> &a
  , const std::vector<T> &b)
{
  if (a.size() != b.size()) throw std::runtime_error("Size mismatch");
  T result = 0.0;
  auto it_b = begin(b);
  for (auto elem_a : a) result += elem_a * (*it_b++);
  return result;
}

/**
 * Calculates zeroth moment of a node based on formula in LBIntro
 * \param node distribution function node containing nine discrete velocity
 *        vectors
 * \return zeroth moment of node
 */
template <typename T>
T GetZerothMoment(const std::vector<T> &node)
{
  T result = 0.0;
  for (auto i : node) result += i;
  return result;
}

/**
 * Calculates first moment of a node based on formula in LBIntro
 * \param node distribution function node containing nine discrete velocity
 *        vectors
 * \return first moment of node
 */
template <typename T>
T GetFirstMoment(const T &node
  , const std::vector<T> &e)
{
  auto nd = e.at(0).size();
  T result(nd, 0);
  for (auto d = 0u; d < nd; ++d) {
    auto it_node = begin(node);
    for (auto dir : e) result[d] += (*it_node++) * dir[d];
  }  // d
  return result;
}

/**
 * Immersed Boundary Method Interpolation Stencil
 * http://lbmworkshop.com/wp-content/uploads/2011/09/2011-08-25_Edmonton_IBM.pdf
 * \param x position in one of the lattice directions
 * \return contribution towards the interpolation function in one of the lattice
 *         directions
 */
template <typename T>
T Phi2(T x)
{
  T phi = 0;
  auto x_abs = fabs(x);
  if (x_abs <= 1) phi = 1 - x_abs;
  return phi;
}

/**
 * Immersed Boundary Method Interpolation Stencil
 * http://lbmworkshop.com/wp-content/uploads/2011/09/2011-08-25_Edmonton_IBM.pdf
 * \param x position in one of the lattice directions
 * \return contribution towards the interpolation function in one of the lattice
 *         directions
 */
template <typename T>
T Phi3(T x)
{
  T phi = 0;
  auto x_abs = fabs(x);
  if (x_abs <= 0.5) {
    phi = (1 + sqrt(1 - 3 * x * x)) / 3;
  }
  else if (x_abs <= 1.5) {
    phi = (5 - 3 * x_abs - sqrt(-2 + 6 * x_abs - 3 * x * x)) / 6;
  }
  return phi;
}

/**
 * Immersed Boundary Method Interpolation Stencil
 * http://lbmworkshop.com/wp-content/uploads/2011/09/2011-08-25_Edmonton_IBM.pdf
 * \param x position in one of the lattice directions
 * \return contribution towards the interpolation function in one of the lattice
 *         directions
 */
template <typename T>
T Phi4(T x)
{
  T phi = 0;
  auto x_abs = fabs(x);
  if (x_abs <= 1) {
    phi = (3 - 2 * x_abs + sqrt(1 + 4 * x_abs - 4 * x * x)) / 8;
  }
  else if (x_abs <= 2) {
    phi = (5 - 2 * x_abs - sqrt(-7 + 12 * x_abs - 4 * x * x)) / 8;
  }
  return phi;
}

/**
 * Approximates Dirac delta-function of IBM using Phi2
 * \param x x-position in lattice
 * \param y y-position in lattice
 * \return approximated Dirac delta-function value
 */

template <typename T>
T Dirac2(T x, T y)
{
  return Phi2(x) * Phi2(y);
}

/**
 * Approximates Dirac delta-function of IBM using Phi3
 * \param x x-position in lattice
 * \param y y-position in lattice
 * \return approximated Dirac delta-function value
 */

template <typename T>
T Dirac3(T x, T y)
{
  return Phi3(x) * Phi3(y);
}

/**
 * Approximates Dirac delta-function of IBM using Phi4
 * \param x x-position in lattice
 * \param y y-position in lattice
 * \return approximated Dirac delta-function value
 */

template <typename T>
T Dirac4(T x, T y)
{
  return Phi4(x) * Phi4(y);
}

#endif  // ALGORITHM_HPP_
