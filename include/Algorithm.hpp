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

template <typename T>
T MagnitudeDifference(const std::vector<T> &a
  , const std::vector<T> &b)
{
  T difference = 0;
  auto it_b = begin(b);
  for (auto i : a) difference += fabs(i - *it_b++) / i;
  return difference;
}

template <typename T, typename U>
bool CheckSteadyState(const T &u_prev
  , const T &u_curr
  , U tolerance)
{
  U diff_sum = 0;
  U sum = 0;
  auto it_u_prev = begin(u_prev);
  for (auto n : u_curr) {
    diff_sum += MagnitudeDifference(n, *it_u_prev++);
  }
  std::cout << diff_sum << std::endl;
  return true;
}
#endif  // ALGORITHM_HPP_
