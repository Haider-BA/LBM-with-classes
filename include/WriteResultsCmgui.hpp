#ifndef WRITERESULTSCMGUI_HPP_
#define WRITERESULTSCMGUI_HPP_

#include <vector>

/**
 * Write the density distribution for display in cmgui. One lbm.exnode file is
 * written out and contains the geometry of the lattice. One lbm<time>.exelem
 * is written out at every time step (note: only one per call to the function)
 * and contains the density/concentration information.
 *
 * \param f The density distribution for the solute concentration
 * \param num_nodes_x The number of lattice points in the x direction
 * \param num_nodes_y The number of lattice points in the y direction
 * \param time The current simulation time (used for filename generation)
 */
void WriteResultsCmgui(
    const std::vector<std::vector<double>> &f
  , int num_nodes_x
  , int num_nodes_y
  , int time);

#endif
