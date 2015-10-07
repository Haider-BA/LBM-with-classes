#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_
#include <vector>
#include "LatticeModel.hpp"
#include "ParticleNode.hpp"

class Particle {
 public:
  /**
   * Constructor: Creates a base particle, containing information such as
   * particle stiffness, number of nodes and center coordinates
   * All units for the particle calculated in physical units instead of lattice
   * units
   * \param stiffness particle stiffness, used to calculate particle force on
   *        fluid
   * \param num_nodes number of boundary nodes in the particle
   * \param center_x particle center's x-coordinate
   * \param center_y particle center's y-coordinate
   * \param mobility toggle to determine whether particle can move in fluid
   * \param lm reference to LatticeModel to provide information on number of
   *        rows, columns, dimensions, discrete directions and lattice velocity
   */
  Particle(double stiffness
    , std::size_t num_nodes
    , double center_x
    , double center_y
    , bool mobility
    , LatticeModel &lm);

  /**
   * Virtual destructor since we are deriving from this class
   */
  virtual ~Particle() = default;

  /**
   * Adds a boundary nodes for the particle
   * \param x x-coordinate
   * \param y y-coordinate
   * \param x_ref reference x-coordinate
   * \param y_ref reference y-coordinate
   * \param u_x node velocity in x-direction
   * \param u_y node velocity in y-direction
   * \param force_x node force in x-direction
   * \param force_y node force in y_direction
   */
  void AddNode(double x
    , double y
    , double x_ref
    , double y_ref
    , double u_x
    , double u_y
    , double force_x
    , double force_y);

  /**
   * Creates a cylinder-shaped particle
   * \param radius radius of particle
   */
  void CreateCylinder(double radius);

  /**
   * Updates the reference positions of the boundary nodes and the center node
   */
  void UpdateReferencePosition();

  /**
   * Enables the rigid particle to move in the fluid
   */
  void ChangeMobility(bool mobility);

  /**
   * Pure virtual method for computing particle forces
   */
  virtual void ComputeForces() = 0;

  /**
   * Returns number of boundary nodes in the particle
   * \return number of nodes in particle
   */
  std::size_t GetNumberOfNodes() const;

  /**
   * Particle node used to store information such as coordinate of particle
   * center
   */
  ParticleNode center;

  /**
   * 1D vector storing boundary nodes of the particle
   */
  std::vector<ParticleNode> nodes;

  /**
   * Particle mobility toggle. Set to TRUE to allow particle to move in fluid.
   * Rigid particles: FALSE by default. Deformable particles: TRUE by default.
   */
  bool is_mobile;

 protected:
  /**
   * Pi, used for calculating areas for cylinder-shaped particles
   */
  double pi_ = 3.14159265;

  /**
   * Particle boundary area
   */
  double area_;

  /**
   * Characteristic lengths of shapes used to recalculate reference node
   * positions. For cylinder, it's radius.
   */
  double char_length_;

  /**
   * Particle stiffness
   */
  double stiffness_;

  /**
   * Number of boundary nodes in particle
   */
  std::size_t number_of_nodes_;

  /**
   * Reference to LatticeModel
   */
  LatticeModel &lm_;
};
#endif  // PARTICLE_HPP_
