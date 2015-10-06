#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "WriteResultsCmgui.hpp"

void WriteResultsCmgui(
    const std::vector<std::vector<double>> &lattice
  , int num_nodes_x
  , int num_nodes_y
  , int time)
{
  static bool first_time = true;
  auto num_nodes = num_nodes_x * num_nodes_y;
  std::vector<double> solute_conc(num_nodes, 0.0);

  // Calculate the solute density
  auto depth = lattice[0].size();

//  for (auto n = 0; n < num_nodes; ++n) {
//    double rho = 0.0;
//    for (auto i = 0u; i < depth; ++i) rho += lattice[n][i];
//    solute_conc[n] = rho;
//  }

//  for (auto n = 0; n < num_nodes; ++n) {
//    double rho = 0.0;
//    for (auto i = 0u; i < depth; ++i) rho += lattice[n][i] * lattice[n][i];
//    solute_conc[n] = sqrt(rho);
//  }

  for (auto n = 0; n < num_nodes; ++n) solute_conc[n] = lattice[n][0];

  // Write out only one 'node' file with the lattice size
  if (first_time) {
    std::ofstream cmgui_node_file;
    cmgui_node_file.open("lbm.exnode");
    cmgui_node_file << " Group name : lbm" << std::endl;
    cmgui_node_file << " #Fields=1" << std::endl;
    cmgui_node_file << " 1) coordinates, coordinate, rectangular cartesian, #Components=2" << std::endl;
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

  // Write out one 'element' file each time step
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
  cmgui_elem_file << " #Fields=2" << std::endl;
  cmgui_elem_file << " 1) coordinates, coordinate, rectangular cartesian, #Components=2" << std::endl;
  cmgui_elem_file << "   x.  l.Lagrange*l.Lagrange, no modify, standard node based." << std::endl;
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
  cmgui_elem_file << "   y.  l.Lagrange*l.Lagrange, no modify, standard node based." << std::endl;
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
  cmgui_elem_file << " 2) solute, field, rectangular cartesian, #Components=1" << std::endl;
  cmgui_elem_file << "   solute.  l.Lagrange*l.Lagrange, no modify, grid based." << std::endl;
  cmgui_elem_file << "   #xi1=" << num_nodes_x -1 << ", #xi2=" << num_nodes_y - 1 << std::endl;
  cmgui_elem_file << " Element:            1 0 0" << std::endl;
  cmgui_elem_file << "   Faces:" << std::endl;
  cmgui_elem_file << "   0 0     1" << std::endl;
  cmgui_elem_file << "   0 0     2" << std::endl;
  cmgui_elem_file << "   0 0     3" << std::endl;
  cmgui_elem_file << "   0 0     4" << std::endl;
  cmgui_elem_file << "   Values:" << std::endl;
  for (auto n : solute_conc) cmgui_elem_file << "  " << n << std::endl;
  cmgui_elem_file << "   Nodes:" << std::endl;
  cmgui_elem_file << "                4            3            2            1" << std::endl;  // BASED ON THE INCORRECT LATTICE ORDER; NEEDS TO BE CHANGED IN OUR CODE, AND THIS SHOULD GO BACK TO 1 2 3 4
  cmgui_elem_file << "   Scale factors:" << std::endl;
  cmgui_elem_file << "       1.0   1.0   1.0   1.0" << std::endl;

  cmgui_elem_file.close();
}
