#include <iostream>
#include "Lattice.hpp"
#include "UnitTest++.h"

SUITE(PrintDemo)
{
static const std::size_t g_num_dimensions = 2;
static const std::size_t g_num_discrete_velocities = 9;
static const std::size_t g_num_rows = 5;
static const std::size_t g_num_cols = 4;
static const double g_dx = 0.0316;
static const double g_dt = 0.001;
static const double g_t_total = 1.0;
static const double g_diffusion_coefficient = 0.2;
static const double g_kinematic_viscosity = 0.2;
static const std::vector<double> g_u0 = {123, 321};
static const std::vector<std::vector<unsigned>> g_src_pos_f = {{0, 0}};
static const std::vector<std::vector<unsigned>> g_src_pos_g = {{0, 0}};
static const std::vector<std::vector<double>> g_src_strength_f = {{1, 1}};
static const std::vector<double> g_src_strength_g = {1};
static const std::vector<std::vector<unsigned>> g_obstacles_pos = {{0, 0}};
static const bool g_is_cd = true;
static const bool g_is_ns = true;
static const bool g_is_not_instant = false;
static const bool g_no_obstacles = false;
static const double g_density_g = 2.0;
static const double g_density_f = 2.0;
static const Lattice g_lattice(g_num_dimensions, g_num_discrete_velocities,
      g_num_rows, g_num_cols, g_dx, g_dt, g_t_total, g_diffusion_coefficient,
      g_kinematic_viscosity, g_density_f, g_density_g, g_u0, g_src_pos_f,
      g_src_pos_g, g_src_strength_f, g_src_strength_g, g_obstacles_pos, g_is_cd,
      g_is_ns, g_is_not_instant, g_no_obstacles);

TEST(PrintDensity)
{
  Lattice lattice(g_lattice);
  lattice.Init(lattice.rho_g, g_density_g);
  std::cout << "rho_g" << std::endl;
  lattice.Print(lattice.rho_g);
}

TEST(PrintEquilibriumDistributionFunction)
{
  Lattice lattice(g_lattice);
  lattice.Init(lattice.rho_g, g_density_g);
  lattice.ComputeEq(lattice.g_eq, lattice.rho_g);
  std::cout << "g_eq" << std::endl;
  lattice.Print(0, lattice.g_eq);
}

TEST(PrintVelocity)
{
  Lattice lattice(g_lattice);
  lattice.Init(lattice.u, g_u0);
  std::cout << "u" << std::endl;
  lattice.Print(1, lattice.u);
}

TEST(PrintBoundary)
{
  Lattice lattice(g_lattice);
  lattice.InitAll();
  lattice.boundary_g = lattice.BoundaryCondition(lattice.g);
  std::cout << "boundary_g" <<std::endl;
  lattice.Print(2, lattice.boundary_g);
}
}  // suite PrintDemo
