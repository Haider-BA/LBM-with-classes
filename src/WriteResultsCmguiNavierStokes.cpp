#include <fstream>
#include <numeric>
#include <vector>
#include "WriteResultsCmguiNavierStokes.hpp"

void WriteResultsCmguiNavierStokes(
    const std::vector<std::vector<double>> &f
  , const std::vector<std::vector<double>> &g
  , const std::vector<bool> &obstacle
  , int num_nodes_x
  , int num_nodes_y
  , int num_nodes
  , int time
  , double density
  , double c
  , double c_s)
{
  static bool first_time = true;
  std::vector<double> temp_u_x(num_nodes);
  std::vector<double> temp_u_y(num_nodes);
  std::vector<double> temp_pressure(num_nodes);
  std::vector<double> temp_solute(num_nodes);
  const double c_s_squared = c_s * c_s;  // square speed of sound

  for (auto y = 0; y < num_nodes_y; ++y) {
    for (auto x = 0; x < num_nodes_x; ++x) {
      const auto n = y*num_nodes_x + x;
      if (obstacle[n]) {  // if obstacle node, nothing is to do
        temp_u_x[n] = 0.0;
        temp_u_y[n] = 0.0;
        // pressure = average pressure
        temp_pressure[n] = density * c_s_squared;
//        temp_solute[n] = g[n][9];
      }
      else {
        // Density, rho, is the sum over the lattice directions, i, of f[i]
        double rho_ns = 0.0;
        for (auto i = 0; i < 9; ++i) rho_ns += f[n][i];
        // x-, and y- velocity components
        temp_u_x[n] = c*(f[n][1] + f[n][5] + f[n][8] -(f[n][3] + f[n][6] +
            f[n][7])) / rho_ns;
        temp_u_y[n] = c*(f[n][2] + f[n][5] + f[n][6] -(f[n][4] + f[n][7] +
            f[n][8])) / rho_ns;
        temp_pressure[n] = rho_ns * c_s_squared;

        double rho_cd = 0.0;
        for (auto i = 0; i < 9; ++i) rho_cd += g[n][i];
        temp_solute[n] = rho_cd;
      }
    }
  }

  if (first_time) {
    std::ofstream cmgui_node_file;
    cmgui_node_file.open("lbm.exnode");
    cmgui_node_file << " Group name : lbm" << std::endl;
    cmgui_node_file << " #Fields=1" << std::endl;
    cmgui_node_file << " 1) coordinates, coordinate, rectangular cartesian, "
        "#Components=2" << std::endl;
    cmgui_node_file << "   x.  Value index= 1, #Derivatives= 0" << std::endl;
    cmgui_node_file << "   y.  Value index= 2, #Derivatives= 0" << std::endl;
    cmgui_node_file << " Node:     1" << std::endl;
    cmgui_node_file << "     0" << std::endl;
    cmgui_node_file << "     0" << std::endl;
    cmgui_node_file << " Node:     2" << std::endl;
    cmgui_node_file << "     " << num_nodes_x - 1 << std::endl;
    cmgui_node_file << "     0" << std::endl;
    cmgui_node_file << " Node:     3" << std::endl;
    cmgui_node_file << "     0" << std::endl;
    cmgui_node_file << "     " << num_nodes_y - 1 << std::endl;
    cmgui_node_file << " Node:     4" << std::endl;
    cmgui_node_file << "     " << num_nodes_x - 1 << std::endl;
    cmgui_node_file << "     " << num_nodes_y - 1 << std::endl;
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
  cmgui_elem_file << "   #xi1=" << num_nodes_x - 1 << ", #xi2=" <<
      num_nodes_y - 1 << std::endl;
  cmgui_elem_file << " 3) u_y, field, rectangular cartesian, #Components=1"
      << std::endl;
  cmgui_elem_file << "   u_y.  l.Lagrange*l.Lagrange, no modify, grid based."
      << std::endl;
  cmgui_elem_file << "   #xi1=" << num_nodes_x - 1 << ", #xi2=" <<
      num_nodes_y - 1 << std::endl;
  cmgui_elem_file << " 4) pressure, field, rectangular cartesian, #Components=1"
      << std::endl;
  cmgui_elem_file << "   pressure.  l.Lagrange*l.Lagrange, no modify, grid "
      "based." << std::endl;
  cmgui_elem_file << "   #xi1=" << num_nodes_x -1 << ", #xi2=" <<
      num_nodes_y - 1 << std::endl;
//  cmgui_elem_file << " 5) obstacle, field, rectangular cartesian, "
//    "#Components=1" << std::endl;
//  cmgui_elem_file << "   obstacle.  l.Lagrange*l.Lagrange, no modify, grid "
//    "based." << std::endl;
//  cmgui_elem_file << "   #xi1=" << num_nodes_x -1 << ", #xi2=" <<
//    num_nodes_y - 1 << std::endl;
//  cmgui_elem_file << " 6) porosity, field, rectangular cartesian, "
//    "#Components=1" << std::endl;
//  cmgui_elem_file << "   porosity.  l.Lagrange*l.Lagrange, no modify, grid "
//    "based." << std::endl;
//  cmgui_elem_file << "   #xi1=" << num_nodes_x -1 << ", #xi2=" << "
//    "num_nodes_y - 1 << std::endl;
  cmgui_elem_file << " 5) solute, field, rectangular cartesian, #Components=1"
      << std::endl;
  cmgui_elem_file << "   solute.  l.Lagrange*l.Lagrange, no modify, grid based."
      << std::endl;
  cmgui_elem_file << "   #xi1=" << num_nodes_x - 1 << ", #xi2=" <<
      num_nodes_y - 1 << std::endl;
  cmgui_elem_file << " Element:            1 0 0" << std::endl;
  cmgui_elem_file << "   Faces:" << std::endl;
  cmgui_elem_file << "   0 0     1" << std::endl;
  cmgui_elem_file << "   0 0     2" << std::endl;
  cmgui_elem_file << "   0 0     3" << std::endl;
  cmgui_elem_file << "   0 0     4" << std::endl;
  cmgui_elem_file << "   Values:" << std::endl;
  for (auto n : temp_u_x) cmgui_elem_file << "  " << n << std::endl;
  for (auto n : temp_u_y) cmgui_elem_file << "  " << n << std::endl;
  for (auto n : temp_pressure) cmgui_elem_file << "  " << n << std::endl;
  for (auto n : obstacle) cmgui_elem_file << "  " << n << std::endl;
//  for (auto n = 0u; n < f.size(); ++n) {
//    cmgui_elem_file << "  " << ((f[n][9] < 1.0) ? "1.0" : "0.0") << std::endl;
//  }
  for (auto n : temp_solute) cmgui_elem_file << "  " << n << std::endl;
  cmgui_elem_file << "   Nodes:" << std::endl;
  cmgui_elem_file << "                1            2            3            4"
      << std::endl;
  cmgui_elem_file << "   Scale factors:" << std::endl;
  cmgui_elem_file << "       1.0   1.0   1.0   1.0" << std::endl;

  cmgui_elem_file.close();
}
