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
 * Checks if velocity (u, v) is in steady state based on the following formula
 * |u_curr - u_prev|/|u_curr| <= tol && |v_curr - v_prev|/|v_curr| <= tol
 * This function knows that the model is 2D. "Dry" boundary nodes such as full-
 * way bounceback nodes need to be excluded from the velocity vectors being
 * passed in.
 * \param u_prev velocity lattice from the previous time step
 * \param u_curr velocity lattice from the current time step
 * \param tolerance absolute tolerance value for steady state checking
 * \return TRUE steady state reached
 *         FALSE steady state not reached
 */
template <typename T, typename U>
bool CheckSteadyState(const T &u_prev
  , const T &u_curr
  , U tolerance)
{
  U u_diff_sum = 0;
  U u_sum = 0;
  U v_diff_sum = 0;
  U v_sum = 0;
  auto nn = u_prev.size();
  for (auto n = 0u; n < nn; ++n) {
    u_diff_sum += fabs(u_curr[n][0] - u_prev[n][0]);
    u_sum += fabs(u_curr[n][0]);
    v_diff_sum += fabs(u_curr[n][1] - u_prev[n][1]);
    v_sum += fabs(u_curr[n][1]);
  }
  return u_diff_sum / u_sum < tolerance && v_diff_sum / v_sum < tolerance;
}
#endif  // ALGORITHM_HPP_
