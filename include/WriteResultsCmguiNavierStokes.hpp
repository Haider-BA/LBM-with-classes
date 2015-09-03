#ifndef WRITERESULTSCMGUINAVIERSTOKES_HPP_
#define WRITERESULTSCMGUINAVIERSTOKES_HPP_

#include <vector>

// NEED TO WRITE SOME DOCUMENTATION HERE
//
// Write out the current fluid velocities in the x and y directions
// along with the pressure field for visualization in cmgui
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
  , double c_s);

#endif
