#include <cmath>  // cos, sin
#include <fstream>
#include <iostream>
#include <vector>
#include "BouncebackNodes.hpp"
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "CollisionNSF.hpp"
#include "ImmersedBoundaryMethod.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeD2Q9.hpp"
#include "ParticleRigid.hpp"
#include "Printing.hpp"
#include "StreamD2Q9.hpp"
#include "StreamPeriodic.hpp"
#include "UnitTest++.h"
#include "WriteResultsCmgui.hpp"
#include "ZouHeNodes.hpp"
#include "ZouHePressureNodes.hpp"

SUITE(SimulationDemo)
{
const static auto g_dx = 0.0316;
const static auto g_dt = 0.001;
const static auto g_cs_sqr = (g_dx / g_dt) * (g_dx / g_dt) / 3.0;
const static auto g_k_visco = 0.2;
const static auto g_rho0_f = 1.0;
const static auto g_d_coeff = 0.2;
const static auto g_rho0_g = 1.0;
const static auto g_time_steps = 1000;
const static auto g_is_instant = true;
const static auto g_is_prestream = true;

TEST(SimulateDiffusion)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_g = {{15, 10}};
  std::vector<double> src_str_g = {50};
  std::vector<double> u0 = {0, 0};
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamPeriodic sp(lm);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g
    , !g_is_instant);
  BouncebackNodes bbcd(lm
    , &cd);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  for (auto t = 0u; t < 501; ++t) {
    g.TakeStep();
    WriteResultsCmgui(g.df, nx, ny, t);
  }
}

TEST(SimulateConvectionDiffusion)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_g = {{15, 10}};
  std::vector<double> src_str_g = {50};
  std::vector<double> u0 = {10, 10};
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamPeriodic sp(lm);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g
    , !g_is_instant);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  for (auto t = 0u; t < 501; ++t) {
    g.TakeStep();
    WriteResultsCmgui(g.df, nx, ny, t);
  }
}

TEST(SimulatePoiseuilleFlow)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {10.0, 0.0});
  std::vector<double> u0 = {0.0, 0.0};
  for (auto n = 0u; n < nx * ny; ++n) src_pos_f.push_back({n % nx, n / nx});
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamPeriodic sp(lm);
  CollisionNSF nsf(lm
    , src_pos_f
    , src_str_f
    , g_k_visco
    , g_rho0_f);
  BouncebackNodes bbnsf(lm
    , &nsf);
  LatticeBoltzmann f(lm
    , nsf
    , sp);
  for (auto x = 0u; x < nx; ++x) {
    bbnsf.AddNode(x, 0);
    bbnsf.AddNode(x, ny - 1);
  }
  f.AddBoundaryNodes(&bbnsf);
  for (auto t = 0u; t < 501; ++t) {
    f.TakeStep();
    WriteResultsCmgui(lm.u, nx, ny, t);
  }
}

TEST(SimulatePoiseuilleFlowHw)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {10.0, 0.0});
  std::vector<double> u0 = {0.0, 0.0};
  for (auto n = 0u; n < nx * ny; ++n) src_pos_f.push_back({n % nx, n / nx});
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamPeriodic sp(lm);
  StreamD2Q9 sd(lm);
  CollisionNSF nsf(lm
    , src_pos_f
    , src_str_f
    , g_k_visco
    , g_rho0_f);
  BouncebackNodes bbnsf(lm
    , &sd);
  LatticeBoltzmann f(lm
    , nsf
    , sd);
  for (auto x = 0u; x < nx; ++x) {
    bbnsf.AddNode(x, 0);
    bbnsf.AddNode(x, ny - 1);
  }
  f.AddBoundaryNodes(&bbnsf);
  for (auto t = 0u; t < 501; ++t) {
    f.TakeStep();
    WriteResultsCmgui(lm.u, nx, ny, t);
  }
}

TEST(SimulateDevelopingPoiseuilleFlow)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<double> u0 = {0.0, 0.0};
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamD2Q9 sd(lm);
  StreamPeriodic sp(lm);
  CollisionNS ns(lm
    , g_k_visco
    , g_rho0_f);
  BouncebackNodes hwbb(lm
    , &sp);
  ZouHeNodes inlet(lm
    , ns);
  ZouHeNodes outlet(lm
    , ns);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  for (auto x = 0u; x < nx; ++x) {
    hwbb.AddNode(x, 0);
    hwbb.AddNode(x, ny - 1);
  }
  for (auto y = 1u; y < ny - 1; ++y) {
    inlet.AddNode(0, y, 0.05, 0.0);
    outlet.AddNode(nx - 1, y, 0.0, 0.0);
  }
  f.AddBoundaryNodes(&inlet);
  f.AddBoundaryNodes(&outlet);
  f.AddBoundaryNodes(&hwbb);
  outlet.ToggleNormalFlow();
  for (auto t = 0u; t < 501; ++t) {
    f.TakeStep();
    WriteResultsCmgui(lm.u, nx, ny, t);
    std::cout << t << std::endl;
  }
}

TEST(SimulateDevelopingPoiseuilleFlowPressureOutlet)
{
  std::size_t ny = 51;
  std::size_t nx = 151;
  auto dx = 0.0316;
  auto dt = dx * dx;
  auto k_visco = 0.005;
  auto rho0_f = 1.0;
  std::vector<double> u0 = {0.0, 0.0};
  LatticeD2Q9 lm(ny
    , nx
    , dx
    , dt
    , u0);
  StreamD2Q9 sd(lm);
  CollisionNS ns(lm
    , k_visco
    , rho0_f);
  BouncebackNodes hwbb(lm
    , &sd);
  BouncebackNodes fwbb(lm
    , &ns);
  ZouHeNodes zhns(lm
    , ns);
  ZouHePressureNodes inlet(lm
    , ns);
  ZouHePressureNodes outlet(lm
    , ns);
  LatticeBoltzmann f(lm
    , ns
    , sd);
  // half-way bounceback on top and bottom wall
  for (auto x = 0u; x < nx; ++x) {
    hwbb.AddNode(x, 0);
    hwbb.AddNode(x, ny - 1);
  }
  auto radius = ny / 4;
  auto x_offset = ny / 2;
  auto y_offset = ny / 2;
  for (auto y = 0u; y < ny; ++y) {
    auto y_pos = abs(y - y_offset);
    for (auto x = 0u; x < nx; ++x) {
      auto x_pos = abs(x - x_offset);
      if (sqrt(y_pos * y_pos + x_pos * x_pos) < radius) fwbb.AddNode(x, y);
    }  // x
  }  // y
  // zou/he velocity inlet, pressure outlet
  for (auto y = 1u; y < ny - 1; ++y) {
    inlet.AddNode(0, y, 1.005);
    outlet.AddNode(nx - 1, y, 1.0);
  }
  f.AddBoundaryNodes(&inlet);
  f.AddBoundaryNodes(&outlet);
  f.AddBoundaryNodes(&hwbb);
  f.AddBoundaryNodes(&fwbb);
  for (auto t = 0u; t < 5001; ++t) {
    std::cout << t << std::endl;
    f.TakeStep();
    if (t % 10 == 0) WriteResultsCmgui(lm.u, nx, ny, t / 10);
  }
}

TEST(SimulateNSCDCoupling)
{
  std::size_t ny = 21;
  std::size_t nx = 31;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f(nx * ny, {50.0, -10.0});
  std::vector<std::vector<std::size_t>> src_pos_g = {{15, 10}};
  std::vector<double> src_str_g = {50};
  std::vector<double> u0 = {-5.0, 0.0};
  for (auto n = 0u; n < nx * ny; ++n) {
    src_pos_f.push_back({n % nx, static_cast<std::size_t>(n / nx)});
  }
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u0);
  StreamPeriodic sp(lm);
  CollisionNSF nsf(lm
    , src_pos_f
    , src_str_f
    , g_k_visco
    , g_rho0_f);
  CollisionCD cd(lm
    , src_pos_g
    , src_str_g
    , g_d_coeff
    , g_rho0_g
    , !g_is_instant);
  BouncebackNodes bbnsf(lm
    , &nsf);
  BouncebackNodes bbcd(lm
    , &cd);
  LatticeBoltzmann f(lm
    , nsf
    , sp);
  LatticeBoltzmann g(lm
    , cd
    , sp);
  for (auto x = 0u; x < nx; ++x) {
    bbnsf.AddNode(x, 0);
    bbnsf.AddNode(x, ny - 1);
    bbcd.AddNode(x, 0);
    bbcd.AddNode(x, ny - 1);
  }
  f.AddBoundaryNodes(&bbnsf);
  g.AddBoundaryNodes(&bbcd);
  for (auto t = 0u; t < 501; ++t) {
    f.TakeStep();
    g.TakeStep();
    WriteResultsCmgui(g.df, nx, ny, t);
  }
}

TEST(SimulateTaylorVortex)
{
  // have to use odd number for sizes
  std::size_t ny = 65;
  std::size_t nx = 65;
  // analytical solution parameters
  std::vector<std::vector<double>> u_lattice_an;
  std::vector<double> rho_lattice_an;
  auto u0_an = 0.001;
  auto k_visco = 0.25;
  auto two_pi = 3.1415926 * 2.0;
  // using one k since it's a square box
  auto k = two_pi / nx;
  for (auto n = 0u; n < nx * ny; ++n) {
    auto x = n % nx;
    auto y = n / nx;
    auto x_an = static_cast<double>(x);
    auto y_an = static_cast<double>(y);
    // analytical formula from "Interpolation methods and the accuracy of
    // lattice-Boltzmann mesh refinement" eq17
    auto u_an = -1.0 * u0_an * cos(k * x_an) * sin(k * y_an);
    auto v_an = u0_an * sin(k * x_an) * cos(k * y_an);
    u_lattice_an.push_back({u_an, v_an});
    auto rho_an = g_rho0_f - 0.25 / g_cs_sqr * u0_an * u0_an *
        (cos(2.0 * k * x_an) + cos(2.0 * k * y_an));
    rho_lattice_an.push_back(rho_an);
  }  // n
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u_lattice_an);
  StreamPeriodic sp(lm);
  CollisionNS ns(lm
    , k_visco
    , rho_lattice_an);
  LatticeBoltzmann f(lm
    , ns
    , sp);
  for (auto t = 0u; t < 501; ++t) {
    f.TakeStep();
    WriteResultsCmgui(lm.u, nx, ny, t);
  }
}

TEST(SimulateTaylorVortexForce)
{
  // have to use odd number for sizes
  std::size_t ny = 65;
  std::size_t nx = 65;
  // analytical solution parameters
  std::vector<std::vector<double>> u_lattice_an;
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f;
  auto u0_an = 0.001;
  auto body_force = sqrt(u0_an);
  auto k_visco = 0.25;
  auto two_pi = 3.1415926 * 2.0;
  // using one k since it's a square box
  auto k = two_pi / nx;
  for (auto n = 0u; n < nx * ny; ++n) {
    auto x = n % nx;
    auto y = n / nx;
    src_pos_f.push_back({x, y});
    auto x_an = static_cast<double>(x);
    auto y_an = static_cast<double>(y);
    // analytical formula from "Interpolation methods and the accuracy of
    // lattice-Boltzmann mesh refinement" eq17
    auto u_an = -1.0 * u0_an * cos(k * x_an) * sin(k * y_an);
    auto v_an = u0_an * sin(k * x_an) * cos(k * y_an);
    u_lattice_an.push_back({u_an, v_an});
    auto f_x = -0.5 / g_cs_sqr * body_force * sin(2.0 * k * x_an);
    auto f_y = -0.5 / g_cs_sqr * body_force * sin(2.0 * k * y_an);
    src_str_f.push_back({f_x, f_y});
  }  // n
  LatticeD2Q9 lm(ny
    , nx
    , g_dx
    , g_dt
    , u_lattice_an);
  StreamPeriodic sp(lm);
  CollisionNSF nsf(lm
    , src_pos_f
    , src_str_f
    , k_visco
    , g_rho0_f);
  LatticeBoltzmann f(lm
    , nsf
    , sp);
  for (auto t = 0u; t < 501; ++t) {
    for (auto n = 0u; n < nx * ny; ++n) {
      auto x = n % nx;
      auto y = n / nx;
      auto x_an = static_cast<double>(x);
      auto y_an = static_cast<double>(y);
      // analytical formula from "Interpolation methods and the accuracy of
      // lattice-Boltzmann mesh refinement" eq17
      auto f_x = -0.5 / g_cs_sqr * body_force * sin(2.0 * k * x_an) * exp(-2.0 *
          k * k * k_visco * t);
      auto f_y = -0.5 / g_cs_sqr * body_force * sin(2.0 * k * y_an) * exp(-2.0 *
          k * k * k_visco * t);
      src_str_f[n] = {f_x, f_y};
    }  // n
    nsf.InitSource(src_pos_f
      , src_str_f);
    f.TakeStep();
    std::cout << t << std::endl;
    WriteResultsCmgui(lm.u, nx, ny, t);
  }
}

TEST(SimulateLidDrivenCavityFlow)
{
  // Reynolds number = velocity * length / viscosity
  std::size_t ny = 256;
  std::size_t nx = 256;
  auto dt = 0.0001;
  auto dx = sqrt(dt);
  std::vector<double> u0 = {0.0, 0.0};
  auto k_visco = 0.080896;
  auto u_lid = 0.316;
  auto v_lid = 0.0;
  LatticeD2Q9 lm(ny
    , nx
    , dx
    , dt
    , u0);
  StreamD2Q9 sd(lm);
  CollisionNS ns(lm
    , k_visco
    , g_rho0_f);
  BouncebackNodes bbns(lm
    , &ns);
  ZouHeNodes zhns(lm
    , ns);
  LatticeBoltzmann f(lm
    , ns
    , sd);
  for (auto y = 0u; y < ny; ++y) {
    bbns.AddNode(0, y);
    bbns.AddNode(nx - 1, y);
  }
  for (auto x = 0u; x < nx; ++x) {
    bbns.AddNode(x, 0);
    zhns.AddNode(x, ny - 1, u_lid, v_lid);
  }
  f.AddBoundaryNodes(&bbns);
  f.AddBoundaryNodes(&zhns);
  for (auto t = 0u; t < 32001; ++t) {
    f.TakeStep();
    if (t % 64 == 0) WriteResultsCmgui(lm.u, nx, ny, t / 64);
    std::cout << t << std::endl;
  }
  std::ofstream myfile;
  myfile.open ("velocities.csv");
  myfile << "u_y,u_x" << std::endl;
  for (auto i = 0u; i < 256; ++i) {
    auto y = 128 * nx + i;
    auto x = i * nx + 128;
    myfile << lm.u[y][1] << "," << lm.u[x][0] << std::endl;
  }
  myfile.close();
}

TEST(SimulateKarmanVortex)
{
  // TODO: try periodic stream
  auto pi = 3.14159265;
  std::size_t nx = 400;
  std::size_t ny = 100;
  auto dt = 0.001;
  auto dx = 0.0316;
  std::vector<double> u0 = {0.0, 0.0};
  std::vector<std::vector<std::size_t>> src_pos_f;
  std::vector<std::vector<double>> src_str_f;
  auto k_visco = 0.03;
  auto u_zh = 0.001;
  auto v_zh = 0.0;
  auto radius = ny / 8 * dx;
  auto stiffness = 1.0 / dx;
  auto center = dx * ny / 2.1;
  std::size_t num_nodes = 2 * pi * radius / 0.3 / dx;
  auto interpolation_stencil = 2;
  LatticeD2Q9 lm(ny
    , nx
    , dx
    , dt
    , u0);
  StreamD2Q9 sd(lm);
  CollisionNSF nsf(lm
    , src_pos_f
    , src_str_f
    , k_visco
    , g_rho0_f);
  BouncebackNodes hwbb(lm
    , &sd);
  ZouHeNodes inlet(lm
    , nsf);
  ZouHeNodes outlet(lm
    , nsf);
  LatticeBoltzmann f(lm
    , nsf
    , sd);
  ParticleRigid cylinder(stiffness
    , num_nodes
    , center
    , center
    , lm);
  cylinder.CreateCylinder(radius);
  ImmersedBoundaryMethod ibm(interpolation_stencil
    , nsf.source
    , lm);
  for (auto x = 0u; x < nx; ++x) {
    hwbb.AddNode(x, 0);
    hwbb.AddNode(x, ny - 1);
  }  // x
  for (auto y = 1u; y < ny - 1; ++y) {
    inlet.AddNode(0, y, u_zh, v_zh);
    outlet.AddNode(nx - 1, y, 0.0, 0.0);
  }  // y
  f.AddBoundaryNodes(&inlet);
  f.AddBoundaryNodes(&outlet);
  f.AddBoundaryNodes(&hwbb);
  ibm.AddParticle(&cylinder);
  outlet.ToggleNormalFlow();
  for (auto node : cylinder.nodes) {
    std::cout << node.coord[0] << " " << node.coord[1] << std::endl;
  }
  auto time = 1001u;
  auto interval = time / 500;
  for (auto t = 0u; t < time; ++t) {
    cylinder.ComputeForces();
    ibm.SpreadForce();
    f.TakeStep();
    ibm.InterpolateFluidVelocity();
    ibm.UpdateParticlePosition();
    if (t % interval == 0) WriteResultsCmgui(lm.u, nx, ny, t / interval);
    std::cout << t << std::endl;
  }
  // cylinder shifts too much
  for (auto node : cylinder.nodes) {
    std::cout << node.coord[0] << " " << node.coord[1] << std::endl;
  }
}
}
