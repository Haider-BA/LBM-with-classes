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
   * Return density lattice for NS
   * \return density lattice for NS
   */
  std::vector<double> GetRhoF() const;

  /**
   * Return density lattice for CD
   * \return density lattice for CD
   */
  std::vector<double> GetRhoG() const;

  std::vector<std::vector<double>> GetVelocity() const;

  void Init(std::vector<bool> &lattice
  , const std::vector<std::vector<std::size_t>> &position);

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
