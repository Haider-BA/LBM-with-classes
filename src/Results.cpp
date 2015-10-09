#include "Results.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "BoundaryNodes.hpp"
#include "LatticeModel.hpp"
#include "LatticeD2Q9.hpp"

Results::Results(LatticeModel &lm)
  : lm_ (lm),
    field_ {3},
    avg_pressure_ {0.0},
    field_names_ {},
    field_nums_ {},
    obstacles_ {}
{
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  obstacles_.assign(nx * ny, false);
  auto result = Results::InitializeCleanFolders();
  if (result != 0) throw std::runtime_error("Error in folder initialization");
}

int Results::InitializeCleanFolders()
{
  // creates output folders if they don't exist
  auto cmgui_folder = system("mkdir -p cmgui_output");
  auto vtk_folder = system("mkdir -p vtk_fluid");
  auto old_cmgui_files = system("rm -f cmgui_output/*");
  auto old_vtk_fluid_files = system("rm -f vtk_fluid/*");
  return cmgui_folder + vtk_folder + old_cmgui_files + old_vtk_fluid_files;
}

void Results::RegisterNS(LatticeBoltzmann *f
  , CollisionModel *ns
  , double initial_density)
{
  const auto c = lm_.GetLatticeSpeed();
  const auto cs_sqr = c * c / 3.0;
  f_ = f;
  ns_ = ns;
  avg_pressure_ = initial_density * cs_sqr;
  ++field_;
  field_nums_.push_back(field_);
  field_names_.push_back("pressure");
  ++field_;
  field_nums_.push_back(field_);
  field_names_.push_back("rho_ns");
}

void Results::RegisterCD(LatticeBoltzmann *g
  , CollisionModel *cd)
{
  g_ = g;
  cd_ = cd;
  ++field_;
  field_nums_.push_back(field_);
  field_names_.push_back("solute");
  ++field_;
  field_nums_.push_back(field_);
  field_names_.push_back("rho_cd");
}

void Results::RegisterObstacles(BoundaryNodes *bn)
{
  for (auto n : bn->position) obstacles_[n] = true;
}

void Results::WriteNode()
{
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  std::ofstream cmgui_node_file;
  cmgui_node_file.open("cmgui_output/lbm.exnode");
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
  cmgui_node_file << "     " << nx - 1 << std::endl;
  cmgui_node_file << "     0" << std::endl;
  cmgui_node_file << " Node:     3" << std::endl;
  cmgui_node_file << "     0" << std::endl;
  cmgui_node_file << "     " << ny - 1 << std::endl;
  cmgui_node_file << " Node:     4" << std::endl;
  cmgui_node_file << "     " << nx - 1 << std::endl;
  cmgui_node_file << "     " << ny - 1 << std::endl;
  cmgui_node_file.close();
}

void Results::WriteResult(int time)
{
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  const auto c = lm_.GetLatticeSpeed();
  const auto cs_sqr = c * c / 3.0;
  std::vector<std::vector<double>> data(2);
  for (auto n : lm_.u) {
    data[0].push_back(n[0]);
    data[1].push_back(n[1]);
  }  // n
  if (ns_) {
    data.push_back(ns_->rho);
    data.push_back(ns_->rho);
    for (auto &i : data[2]) i *= cs_sqr;
  }
  if (cd_) {
    data.push_back(std::vector<double>());
    const auto m = data.size() - 1;
    for (auto n : g_->df) {
      auto temp = 0.0;
      for (auto i : n) temp += i;
      data[m].push_back(temp);
    }  // n
    data.push_back(cd_->rho);
  }
  for (auto n = 0u; n < nx * ny; ++n) {
    if (obstacles_[n]) {
      for (auto m = 0u; m < data.size(); ++m) {
        data[m][n] = (ns_ && m == 2) ? avg_pressure_ : 0.0;
      }  // m
    }
  }  // n
  std::ofstream cmgui_elem_file;
  cmgui_elem_file.open("cmgui_output/lbm" + std::to_string(time) + ".exelem");
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
  cmgui_elem_file << " #Fields=" << field_ << std::endl;
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
  cmgui_elem_file << " 2) u_x, field, rectangular cartesian, #Components=1" <<
      std::endl;
  cmgui_elem_file << "   u_x.  l.Lagrange*l.Lagrange, no modify, grid based." <<
      std::endl;
  cmgui_elem_file << "   #xi1=" << nx - 1 << ", #xi2=" << ny - 1 << std::endl;
  cmgui_elem_file << " 3) u_y, field, rectangular cartesian, #Components=1" <<
      std::endl;
  cmgui_elem_file << "   u_y.  l.Lagrange*l.Lagrange, no modify, grid based." <<
      std::endl;
  cmgui_elem_file << "   #xi1=" << nx - 1 << ", #xi2=" << ny - 1 << std::endl;
  for (auto i = 0u; i < field_nums_.size(); ++i) {
    cmgui_elem_file << " " << field_nums_[i] << ") " + field_names_[i] + ", "
        "field, rectangular cartesian, #Components=1" << std::endl;
    cmgui_elem_file << "   " + field_names_[i] + ".  l.Lagrange*l.Lagrange, no "
        "modify, grid based." << std::endl;
    cmgui_elem_file << "   #xi1=" << nx - 1 << ", #xi2=" << ny - 1 << std::endl;
  }
  cmgui_elem_file << " Element:            1 0 0" << std::endl;
  cmgui_elem_file << "   Faces:" << std::endl;
  cmgui_elem_file << "   0 0     1" << std::endl;
  cmgui_elem_file << "   0 0     2" << std::endl;
  cmgui_elem_file << "   0 0     3" << std::endl;
  cmgui_elem_file << "   0 0     4" << std::endl;
  cmgui_elem_file << "   Values:" << std::endl;
  // write secondary data, ns: pressure, rho_ns, cd: solute, rho_cd
  for (auto i = 0; i < field_ - 1; ++i) {
    for (auto n : data[i]) cmgui_elem_file << "  " << n << std::endl;
  }  // i
  cmgui_elem_file << "   Nodes:" << std::endl;
  cmgui_elem_file << "                1            2            3            4"
      << std::endl;
  cmgui_elem_file << "   Scale factors:" << std::endl;
  cmgui_elem_file << "       1.0   1.0   1.0   1.0" << std::endl;
  cmgui_elem_file.close();
}

void Results::WriteResultVTK(int time)
{
  if (!ns_) throw std::runtime_error("Not NS");
  const auto nx = lm_.GetNumberOfColumns();
  const auto ny = lm_.GetNumberOfRows();
  std::ofstream vtk_file;
  vtk_file.open("vtk_fluid/fluid_t" + std::to_string(time) + ".vtk");
  vtk_file << "# vtk DataFile Version 3.0" << std::endl;
  vtk_file << "fluid_state" << std::endl;
  vtk_file << "ASCII" << std::endl;
  vtk_file << "DATASET RECTILINEAR_GRID" << std::endl;
  vtk_file << "DIMENSIONS " << nx << " " << ny << " 1" << std::endl;

  // Write x, y, z coordinates. z set to be 1 since it's 2D
  vtk_file << "X_COORDINATES " << nx << " float" << std::endl;
  for (auto x = 0u; x < nx; ++x) vtk_file << x << " ";
  vtk_file << std::endl;
  vtk_file << "Y_COORDINATES " << ny << " float" << std::endl;
  for (auto y = 0u; y < ny; ++y) vtk_file << y << " ";
  vtk_file << std::endl;
  vtk_file << "Z_COORDINATES " << 1 << " float" << std::endl;
  vtk_file << 0 << std::endl;
  vtk_file << "POINT_DATA " << nx * ny << std::endl;

  // Write density difference
  vtk_file << "SCALARS density_difference float 1" << std::endl;
  vtk_file << "LOOKUP_TABLE default" << std::endl;
  for (auto density : ns_->rho) vtk_file << density - 1.0 << std::endl;

  // Write velocity as vectors
  vtk_file << "VECTORS velocity_vector float" << std::endl;
  for (auto v : lm_.u) vtk_file << v[0] << " " << v[1] << " 0" << std::endl;
  vtk_file.close();
}
