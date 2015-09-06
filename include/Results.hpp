#ifndef RESULTS_HPP_
#define RESULTS_HPP_

#include <string>
#include <vector>
#include "BoundaryNodes.hpp"
#include "CollisionModel.hpp"
#include "LatticeBoltzmann.hpp"
#include "LatticeModel.hpp"
#include "LatticeD2Q9.hpp"

class Results {
 public:
  /**
   * Constructor: Creates results class with reference to LatticeModel for
   * information on number of rows, columns, space step, time step and lattice
   * velocity
   * \param lm reference to LatticeModel
   */
  Results(LatticeModel &lm);

  /**
   * Override copy constructor due to -Weffc++ warnings
   */
  Results(const Results&) = default;

  /**
   * Override copy assignment due to -Weffc++ warnings
   */
  Results& operator= (const Results&) = default;

  /**
   * Destructor
   */
  ~Results() = default;

  /**
   * Registers information about Navier-Stokes equation to results
   * \param f pointer to distribution functions for Navier-Stokes equation
   * \param ns pointer to collision model for Navier-Stokes equation, contains
   *        information on density (pressure)
   * \param initial_density initial density of f, for use when there are
   *        obstacles
   */
  void RegisterNS(LatticeBoltzmann *f
    , CollisionModel *ns
    , double initial_density);

  /**
   * Registers information about Convection-diffusion equation to results
   * \param g pointer to distribution functions for Convection-diffusion
   *        equation
   * \param cd pointer to collision model for Convection-diffusion equation,
   *        contains information on density (concentration of solute)
   */
  void RegisterCD(LatticeBoltzmann *g
    , CollisionModel *cd);

  /**
   * Registers obstacles positions in the lattice
   * \param bn pointer to obstacles
   */
  void RegisterObstacles(BoundaryNodes *bn);

  /**
   * Writes .exnode file
   */
  void WriteNode();

  /**
   * Writes results at a particular time point. Currently writes: coordinates,
   * velocity in x- and y- direction, pressure and density of NS (if NS is
   * present), solute concentration and density of CD (if CD is present)
   * \param time time point
   * \return void
   *
   */
  void WriteResult(int time);

 private:
  /**
   * Reference to LatticeModel
   */
  LatticeModel &lm_;

  /**
   * Total number of fields in output
   */
  int field_;

  /**
   * Average pressure (in NS)
   */
  double avg_pressure_;

  /**
   * Vector containing field names, based on what information are registered in
   * Results
   */
  std::vector<std::string> field_names_;

  /**
   * Vector containing field number, based on what information are registered in
   * Results
   */
  std::vector<int> field_nums_;

  /**
   * Vector containing Boolean values to indicate presence of obstacles at a
   * node
   */
  std::vector<bool> obstacles_;

  /**
   * Pointer to LBM class for Navier-Stokes equation
   */
  LatticeBoltzmann *f_ = nullptr;

  /**
   * Pointer to LBM class for Convection-Diffusion equation
   */
  LatticeBoltzmann *g_ = nullptr;

  /**
   * Pointer to CollisionModel class for Navier-Stokes equation
   */
  CollisionModel *ns_ = nullptr;

  /**
   * Pointer to CollisionModel class for Convection-Diffusion equation
   */
  CollisionModel *cd_ = nullptr;
};

#endif  // RESULTS_HPP_
