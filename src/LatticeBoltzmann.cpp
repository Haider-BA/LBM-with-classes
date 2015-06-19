#include "LatticeBoltzmann.hpp"
#include <iostream>
#include <stdexcept>  // runtime_error

LatticeBoltzmann::LatticeBoltzmann()
{
  throw std::runtime_error("Declaring uninitialized LBM is not allowed");
}

LatticeBoltzmann::LatticeBoltzmann(std::size_t num_dims
  , std::size_t num_dirs
  , std::size_t num_rows
  , std::size_t num_cols
  , double dx
  , double dt
  , double t_total
  , double diffusion_coefficient
  , double kinematic_viscosity
  , double initial_density_f
  , double initial_density_g
  , const std::vector<double> &u0
  , const std::vector<std::vector<std::size_t>> &src_position_f
  , const std::vector<std::vector<double>> &src_strength_f
  , const std::vector<std::vector<std::size_t>> &src_position_g
  , const std::vector<double> &src_strength_g
  , const std::vector<std::vector<std::size_t>> &obstacles_position
  , bool is_ns
  , bool is_cd
  , bool is_instant
  , bool has_obstacles)
  : number_of_dimensions_ {num_dims},
    number_of_directions_ {num_dirs},
    number_of_rows_ {num_rows},
    number_of_columns_ {num_cols},
    dx_ {dx},
    dt_ {dt},
    total_time_ {t_total},
    is_ns_ {is_ns},
    is_cd_ {is_cd},
    is_instant_ {is_instant},
    has_obstacles_ {has_obstacles}
{
  if (CheckParameters()) {
    throw std::runtime_error("Zero value in input parameters");
  }
  else {
    // Initializing all variables with zero values
    std::size_t lattice_size = (num_rows) * (num_cols);
    std::vector<double> length_d(num_dims, 0.0);
    std::vector<double> length_q(num_dirs, 0.0);
    // Initializing variables with initial values
    c_ = dx_ / dt_;
    cs_sqr_ = c_ * c_ / 3.0;
    u.assign(lattice_size, u0);
    if (is_ns_) {
      f.assign(lattice_size, length_q);
      f_eq.assign(lattice_size, length_q);
      rho_f.assign(lattice_size, initial_density_f);
      boundary_f.assign(2 * num_rows + 2 * num_cols + 4, length_q);
      src_f.assign(lattice_size, length_d);
      // tau_ns_ formula from "Discrete lattice effects on the forcing term in
      // the lattice Boltzmann method" Guo2002
      tau_ns_ = 0.5 + kinematic_viscosity / cs_sqr_ / dt_;
    }
    if (is_cd_) {
      g.assign(lattice_size, length_q);
      g_eq.assign(lattice_size, length_q);
      rho_g.assign(lattice_size, initial_density_g);
      boundary_g.assign(2 * num_rows + 2 * num_cols + 4, length_q);
      src_g.assign(lattice_size, 0.0);
      // tau_cd_ formula from "A new scheme for source term in LBGK model for
      // convection diffusion equation"
      tau_cd_ = 0.5 + diffusion_coefficient / cs_sqr_ / dt_;
    }
    if (has_obstacles) {
      obstacles.assign(lattice_size, false);
      LatticeBoltzmann::Init(obstacles, obstacles_position);
    }
  }
}

std::size_t LatticeBoltzmann::GetNumberOfDimensions() const
{
  return number_of_dimensions_;
}

std::size_t LatticeBoltzmann::GetNumberOfDirections() const
{
  return number_of_directions_;
}

std::size_t LatticeBoltzmann::GetNumberOfRows() const
{
  return number_of_rows_;
}

std::size_t LatticeBoltzmann::GetNumberOfColumns() const
{
  return number_of_columns_;
}

void LatticeBoltzmann::Init(std::vector<bool> &lattice
  , const std::vector<std::vector<std::size_t>> &position)
{
  auto nx = GetNumberOfColumns();
  auto ny = GetNumberOfRows();
  auto nd = GetNumberOfDimensions();
  if(position[0].size() != nd)
      throw std::runtime_error("Insufficient position information");
  for (auto pos : position) {
    if (pos[0] > nx - 1) throw std::runtime_error("x value out of range");
    if (pos[1] > ny - 1) throw std::runtime_error("y value out of range");
    lattice[pos[1] * nx + pos[0]] = true;
  }  // pos
}

bool LatticeBoltzmann::CheckParameters()
{
  return number_of_dimensions_ == 0 || number_of_directions_ == 0 ||
      number_of_rows_ == 0 || number_of_columns_ == 0 || dx_ < 1e-20 ||
      dt_ < 1e-20 || total_time_ < 1e-20 || !(is_ns_ || is_cd_);
}
