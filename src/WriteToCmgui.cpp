#include "WriteToCmgui.hpp"
#include <fstream>
#include <numeric>
#include <vector>

// TODO: change coordinate system

void WriteToCmgui(const std::vector<std::vector<double>> &f
  , const std::vector<std::vector<double>> &g
  , const std::vector<std::vector<double>> &u
  , std::size_t nx
  , std::size_t ny
  , int time
  , double rho0
  , double c
  , double cs_sqr)
{
  enum Directions {
    E = 1,
    N,
    W,
    S,
    NE,
    NW,
    SW,
    SE
  };
  static bool first_time = true;
  auto num_nodes = nx * ny;
  std::vector<double> u_x(num_nodes);
  std::vector<double> u_y(num_nodes);
  std::vector<double> pressure(num_nodes);
  std::vector<double> solute(num_nodes);

  for (auto y = 0u; y < ny; ++y) {
    for (auto x = 0u; x < nx; ++x) {
      const auto n = y * nx + x;
      // add code to detect boundary nodes
      auto rho_ns = 0.0;
      // density is sum of all discrete velocities per node
      for (auto i : f[n]) rho_ns += i;
      u_x[n] = c * (f[n][E] + f[n][NE] + f[n][SE] - f[n][W] - f[n][NW] -
          f[n][SW]) / rho_ns;
      u_y[n] = c * (f[n][N] + f[n][NE] + f[n][NW] - f[n][S] - f[n][SW] -
          f[n][SE]) / rho_ns;
      pressure[n] = rho_ns * cs_sqr;

      auto rho_cd = 0.0;
      for (auto i : g[n]) rho_cd += i;
      solute[n] = rho_cd;
    }
  }

  if (first_time) {
    std::ofstream cmgui_node_file;
    cmgui_node_file.open("lbm.exnode");
    cmgui_node_file << " Group name : lbm" << std::endl;
    cmgui_node_file << " #Fields=1" << std::endl;
    cmgui_node_file << " 1) coordinates, coordinate, rectangular cartesian, #Co"
        "mponents=2" << std::endl;
    cmgui_node_file << "   x.  Value index= 1, #Derivatives= 0" << std::endl;
    cmgui_node_file << "   y.  Value index= 2, #Derivatives= 0" << std::endl;
    cmgui_node_file << " Node:     1" << std::endl;
    cmgui_node_file << "     0" << std::endl;
    cmgui_node_file << "     0" << std::endl;
    cmgui_node_file << " Node:     2" << std::endl;
    cmgui_node_file << "     " << nx - 1 << std::endl;
    cmgui_node_file << "     0" << std::endl;
    cmgui_node_file << " Node:     3" << std::endl;
    cmgui_node_file << "     0" << std::endl;
    cmgui_node_file << "     " << ny - 1 << std::endl;
    cmgui_node_file << " Node:     4" << std::endl;
    cmgui_node_file << "     " << nx - 1 << std::endl;
    cmgui_node_file << "     " << ny - 1 << std::endl;
    cmgui_node_file.close();
    first_time = false;
  }

  std::ofstream cmgui_elem_file;
  cmgui_elem_file.open("lbm"+std::to_string(time)+".exelem");
  cmgui_elem_file << " Group name: lbm" << std::endl;
  cmgui_elem_file << " Shape.  Dimension=1" << std::endl;
  cmgui_elem_file << " Element: 0 0     1" << std::endl;
  cmgui_elem_file << " Element: 0 0     2" << std::endl;
  cmgui_elem_file << " Element: 0 0     3" << std::endl;
  cmgui_elem_file << " Element: 0 0     4" << std::endl;
  cmgui_elem_file << " Shape.  Dimension=2" << std::endl;
  cmgui_elem_file << " #Scale factor sets= 1" << std::endl;
  cmgui_elem_file << "   l.Lagrange*l.Lagrange, #Scale factors= 4" << std::endl;
  cmgui_elem_file << " #Nodes=           4" << std::endl;
  cmgui_elem_file << " #Fields=5" << std::endl;
  cmgui_elem_file << " 1) coordinates, coordinate, rectangular cartesian, "
      "#Components=2" << std::endl;
  cmgui_elem_file << "   x.  l.Lagrange*l.Lagrange, no modify, standard node "
      "based." << std::endl;
  cmgui_elem_file << "     #Nodes= 4" << std::endl;
  cmgui_elem_file << "      1.  #Values=1" << std::endl;
  cmgui_elem_file << "       Value indices:     1" << std::endl;
  cmgui_elem_file << "       Scale factor indices:   1" << std::endl;
  cmgui_elem_file << "      2.  #Values=1" << std::endl;
  cmgui_elem_file << "       Value indices:     1" << std::endl;
  cmgui_elem_file << "       Scale factor indices:   2" << std::endl;
  cmgui_elem_file << "      3.  #Values=1" << std::endl;
  cmgui_elem_file << "       Value indices:     1" << std::endl;
  cmgui_elem_file << "       Scale factor indices:   3" << std::endl;
  cmgui_elem_file << "      4.  #Values=1" << std::endl;
  cmgui_elem_file << "       Value indices:     1" << std::endl;
  cmgui_elem_file << "       Scale factor indices:   4" << std::endl;
  cmgui_elem_file << "   y.  l.Lagrange*l.Lagrange, no modify, standard node "
      "based." << std::endl;
  cmgui_elem_file << "     #Nodes= 4" << std::endl;
  cmgui_elem_file << "      1.  #Values=1" << std::endl;
  cmgui_elem_file << "       Value indices:     1" << std::endl;
  cmgui_elem_file << "       Scale factor indices:   1" << std::endl;
  cmgui_elem_file << "      2.  #Values=1" << std::endl;
  cmgui_elem_file << "       Value indices:     1" << std::endl;
  cmgui_elem_file << "       Scale factor indices:   2" << std::endl;
  cmgui_elem_file << "      3.  #Values=1" << std::endl;
  cmgui_elem_file << "       Value indices:     1" << std::endl;
  cmgui_elem_file << "       Scale factor indices:   3" << std::endl;
  cmgui_elem_file << "      4.  #Values=1" << std::endl;
  cmgui_elem_file << "       Value indices:     1" << std::endl;
  cmgui_elem_file << "       Scale factor indices:   4" << std::endl;
  cmgui_elem_file << " 2) u_x, field, rectangular cartesian, #Components=1"
      << std::endl;
  cmgui_elem_file << "   u_x.  l.Lagrange*l.Lagrange, no modify, grid based."
      << std::endl;
  cmgui_elem_file << "   #xi1=" << nx - 1 << ", #xi2=" << ny - 1 << std::endl;
  cmgui_elem_file << " 3) u_y, field, rectangular cartesian, #Components=1"
      << std::endl;
  cmgui_elem_file << "   u_y.  l.Lagrange*l.Lagrange, no modify, grid based."
      << std::endl;
  cmgui_elem_file << "   #xi1=" << nx - 1 << ", #xi2=" << ny - 1 << std::endl;
  cmgui_elem_file << " 4) pressure, field, rectangular cartesian, #Components=1"
      << std::endl;
  cmgui_elem_file << "   pressure.  l.Lagrange*l.Lagrange, no modify, grid "
      "based." << std::endl;
  cmgui_elem_file << "   #xi1=" << nx -1 << ", #xi2=" << ny - 1 << std::endl;
  cmgui_elem_file << " 5) solute, field, rectangular cartesian, #Components=1"
      << std::endl;
  cmgui_elem_file << "   solute.  l.Lagrange*l.Lagrange, no modify, grid based."
      << std::endl;
  cmgui_elem_file << "   #xi1=" << nx -1 << ", #xi2=" << ny - 1 << std::endl;
  cmgui_elem_file << " Element:            1 0 0" << std::endl;
  cmgui_elem_file << "   Faces:" << std::endl;
  cmgui_elem_file << "   0 0     1" << std::endl;
  cmgui_elem_file << "   0 0     2" << std::endl;
  cmgui_elem_file << "   0 0     3" << std::endl;
  cmgui_elem_file << "   0 0     4" << std::endl;
  cmgui_elem_file << "   Values:" << std::endl;
  for (auto n : u_x) cmgui_elem_file << "  " << n << std::endl;
  for (auto n : u_y) cmgui_elem_file << "  " << n << std::endl;
  for (auto n : pressure) cmgui_elem_file << "  " << n << std::endl;
  for (auto n : solute) cmgui_elem_file << "  " << n << std::endl;
  cmgui_elem_file << "   Nodes:" << std::endl;
  cmgui_elem_file << "                4            3            2            1"
      << std::endl;
  cmgui_elem_file << "   Scale factors:" << std::endl;
  cmgui_elem_file << "       1.0   1.0   1.0   1.0" << std::endl;

  cmgui_elem_file.close();
}
