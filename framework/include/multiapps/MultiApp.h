#ifndef MULTIAPP_H
#define MULTIAPP_H

#include "MooseApp.h"

// libMesh includes
#include "libmesh/mesh_tools.h"

class MultiApp;

template<>
InputParameters validParams<MultiApp>();

/**
 * A MultiApp represents one or more MOOSE applications that are running simultaneously.
 * These other MOOSE apps generally represent some "sub-solve" or "embedded-solves"
 * of the overall nonlinear solve.
 */
class MultiApp : public MooseObject
{
public:
  MultiApp(const std::string & name, InputParameters parameters);

  virtual ~MultiApp();

  /**
   * Re-solve all of the Apps.
   */
  virtual void solveStep() = 0;

  /**
   * @param app The global app number to get the Executioner for
   * @return The Executioner associated with that App.
   */
  virtual Executioner * getExecutioner(unsigned int app);

  /**
   * Get the BoundingBox for the mesh associated with app
   * @param app The global app number you want to get the bounding box for
   */
  virtual MeshTools::BoundingBox getBoundingBox(unsigned int app);

  /**
   * @return When this MultiApp will be executed.
   */
  int executeOn() { return _execute_on; }

  /**
   * Get the FEProblem this MultiApp is part of.
   */
  FEProblem * problem() { return _fe_problem; }

  /**
   * Get the FEProblem for the global app is part of.
   * @param app The global app number
   */
  FEProblem * appProblem(unsigned int app);

  /**
   * Get a UserObject base for a specific global app
   * @param app The global app number you want to get a UserObject from.
   */
  const UserObject & appUserObjectBase(unsigned int app, const std::string & name);

  /**
   * @return Number of Apps on local processor.
   */
  unsigned int numLocalApps() { return _my_num_apps; }

  /**
   * @return The global number of the first app on the local processor.
   */
  unsigned int firstLocalApp() { return _first_local_app; }

  /**
   * The physical position of a global App number
   * @param app The global app number you want the position for.
   * @return the position
   */
  Point position(unsigned int app) { return _positions[app]; }

protected:
  /**
   * Create an MPI communicator suitable for all of these apps.
   */
  void buildComm();

  /**
   * Swap the libMesh MPI communicator out for ours.
   */
  MPI_Comm swapLibMeshComm(MPI_Comm new_comm);

  /**
   * Map a global App number to the local number.
   * Note: This will error if given a global number that doesn't map to a local number.
   *
   * @param global_app The global app number.
   * @return The local app number.
   */
  unsigned int globalAppToLocal(unsigned int global_app);

  /// The FEProblem this MultiApp is part of
  FEProblem * _fe_problem;

  /// The type of application to build
  std::string _app_type;

  /// The positions as they came in from the input file
  std::vector<Real> _positions_vec;

  /// The positions of all of the apps
  std::vector<Point> _positions;

  /// The input file for each app's simulation
  std::vector<std::string> _input_files;

  /// The output file basename for each multiapp
  std::string _output_base;

  /// The total number of apps to simulate
  unsigned int _total_num_apps;

  /// The number of apps this object is involved in simulating
  unsigned int _my_num_apps;

  /// The number of the first app on this processor
  unsigned int _first_local_app;

  /// The comm that was passed to us specifying our pool of processors
  MPI_Comm _orig_comm;

  /// The MPI communicator this object is going to use.
  MPI_Comm _my_comm;

  /// The number of processors in the origal comm
  unsigned int _orig_num_procs;

  /// The mpi "rank" of this processor in the original communicator
  unsigned int _orig_pid;

  /// Pointers to each of the Apps
  std::vector<MooseApp *> _apps;

  /// When this MultiApp will be executed
  MooseEnum _execute_on;
};

#endif // MULTIAPP_H
