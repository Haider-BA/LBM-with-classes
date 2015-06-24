#ifndef LATTICEBOLTZMANN_HPP_
#define LATTICEBOLTZMANN_HPP_
#include <vector>
#include "CollisionCD.hpp"
#include "CollisionNS.hpp"
#include "LatticeModel.hpp"

class LatticeBoltzmann {
 public:
  /**
   * Constructor: (default) Override default constructor to throw exception to
   * to forbid declaring uninitialized
   * Does not work because of uninitialized reference member error (lm_)
   * Can try calling intializing lm with 0's
   */
//  LatticeBoltzmann();

  /**
   * Constructor: Creates lattice
   * \param
   * \param
   * \return
   *
   */
  LatticeBoltzmann(double t_total
    , const std::vector<std::vector<std::size_t>> &obstacles_pos
    , bool is_ns
    , bool is_cd
    , bool is_instant
    , bool has_obstacles
    , LatticeModel &lm
    , CollisionNS &ns
    , CollisionCD &cd);

  ~LatticeBoltzmann() = default;

  LatticeBoltzmann(const LatticeBoltzmann&) = default;

   /**
   * Get the number of dimensions of the lattice. 2 for 2D and 3 for 3D.
   *
   * \return number of dimensions of the lattice
   */
  std::size_t GetNumberOfDimensions() const;

  /**
   * Get the number of discrete velocities of the lattice, specified by the
   * model used. 9 for Q9.
   *
   * \return number of discrete velocities of the lattice
   */
  std::size_t GetNumberOfDirections() const;

  /**
   * Get the number of rows of the lattice, will be 2 smaller than the
   * actual number of rows created due to additional boundary layer.
   *
   * \return number of rows of the lattice
   */
  std::size_t GetNumberOfRows() const;

  /**
   * Get the number of columns of the lattice, will be 2 smaller than
   * the actual number of columns created due to additional boundary layer.
   *
   * \return number of columns of the lattice
   */
  std::size_t GetNumberOfColumns() const;

  /**
   * Initializes obstacle lattice
   * \param lattice obstacles at each node stored row-wise
   * \param position 2D vector containing position information of the obstacles
   */
  void Init(std::vector<bool> &lattice
    , const std::vector<std::vector<std::size_t>> &position);

  /**
   * Reflects distribution functions on obstacle nodes, checks if node is not at
   * edges of the reflecting direction and if node in the reflecting direction
   * is not an obstacle node before reflecting
   * \param lattice 2D vector containing distribution functions
   */
  void Obstacles(std::vector<std::vector<double>> &lattice);

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  std::vector<std::vector<double>> BoundaryCondition(
      const std::vector<std::vector<double>> &lattice);

  /** \brief
   *
   * \param
   * \param
   * \return
   *
   */
  std::vector<std::vector<double>> Stream(
      const std::vector<std::vector<double>> &lattice
    , const std::vector<std::vector<double>> &boundary);

  /** \brief
   *
   * \return void
   *
   */
  void TakeStep();

  /** \brief
   *
   * \return void
   *
   */
  void RunSim();

  /** \brief
   *
   * \param lattice std::vector<std::vector<double>>&
   * \return void
   *
   */
  void RunSim(std::vector<std::vector<double>> &lattice);

  std::vector<double> Flip(const std::vector<double> &lattice);
  std::vector<bool> Flip(const std::vector<bool> &lattice);
  std::vector<std::vector<double>> Flip(
      const std::vector<std::vector<double>> &lattice);
  void Print(const std::vector<double> &lattice);
  void Print(const std::vector<bool> &lattice);
  void Print(int which_to_print
  , const std::vector<std::vector<double>> &lattice);

  /**
   * NS distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> f;

  /**
   * Boundary nodes for NS lattice
   */
  std::vector<std::vector<double>> boundary_f;

  /**
   * CDE distribution function stored row-wise in a 2D vector.
   */
  std::vector<std::vector<double>> g;

  /**
   * Boundary nodes for CDE lattice
   */
  std::vector<std::vector<double>> boundary_g;

  /**
   * Lattice containing obstacles
   */
  std::vector<bool> obstacles;

 private:
  // 6  2  5  ^
  //  \ | /   |
  // 3--0--1  |
  //  / | \   |
  // 7  4  8  +------->
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
  /**
   * Checks input parameters to ensure there's not invalid values
   * \return true if there is invalid: 0 in any values, both is_ns_ and is_cd_
   *         are false
   *         false if all input values are valid
   */
  bool CheckParameters();
  // input parameters
  double total_time_;
  bool is_ns_;
  bool is_cd_;
  bool is_instant_;
  bool has_obstacles_;
  // LatticeModel to take care of dims, dirs, rows, cols and discrete e vectors
  // by reference, similar to by pointer
  // https://stackoverflow.com/questions/9285627/is-it-possible-to-pass-derived-
  // classes-by-reference-to-a-function-taking-base-cl
  LatticeModel &lm_;
  // Collision models to take care of eq, rho, source
  CollisionNS &ns_;
  CollisionCD &cd_;
};
#endif // LATTICEBOLTZMANN_HPP_
