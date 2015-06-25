#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_
#include <stdexcept>
#include <vector>
#include <iostream>
#include <typeinfo>

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
std::vector<T> GetFirstMoment(const std::vector<T> &node
  , const std::vector<std::vector<T>> &e)
{
  auto nd = e.at(0).size();
  std::vector<T> result(nd, 0);
  for (auto d = 0u; d < nd; ++d) {
    auto it_node = begin(node);
    for (auto dir : e) result[d] += (*it_node++) * dir[d];
  }  // d
  return result;
}

#endif  // ALGORITHM_HPP_
