#ifndef AUXKERNEL_H
#define AUXKERNEL_H

#include "MooseObject.h"
#include "Coupleable.h"
#include "MaterialPropertyInterface.h"
#include "FunctionInterface.h"
#include "TransientInterface.h"
#include "PostprocessorInterface.h"
#include "GeometricSearchInterface.h"
#include "MooseMesh.h"

//forward declarations
class Problem;
class SubProblemInterface;
class AuxKernel;
class AuxiliarySystem;

template<>
InputParameters validParams<AuxKernel>();

/** 
 * AuxKernels compute values at nodes.
 */
class AuxKernel :
  public MooseObject,
  public Coupleable,
  public FunctionInterface,
  public TransientInterface,
  public MaterialPropertyInterface,
  public PostprocessorInterface,
  protected GeometricSearchInterface
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  AuxKernel(const std::string & name, InputParameters parameters);

  virtual ~AuxKernel() {}

  void compute(NumericVector<Number> & sln);

  MooseVariable & variable() { return _var; }

  virtual void setup() { }

  bool isNodal();

  bool ts() { return _ts; }

protected:
  virtual Real computeValue() = 0;

  Problem & _problem;
  SubProblemInterface & _subproblem;
  AuxiliarySystem & _aux_sys;
  THREAD_ID _tid;
  MooseVariable & _var;
  MooseMesh & _mesh;
  int _dim;

  const std::vector< Point > & _q_point;
  QBase * & _qrule;
  const std::vector<Real> & _JxW;

  VariableValue & _u;                                   /// Holds the solution at current quadrature points
  VariableValue & _u_old;                               /// Holds the previous solution at the current quadrature point.
  VariableValue & _u_older;                             /// Holds the t-2 solution at the current quadrature point.

  const Elem * & _current_elem;
  const Node * & _current_node;

  Real & _current_volume;

  /**
   * true if the kernel nodal, false if it is elemental
   */
  bool _nodal;

  unsigned int _qp;

  bool _ts;

  // Single Instance Variables
  Real & _real_zero;
  Array<Real> & _zero;
  Array<RealGradient> & _grad_zero;
  Array<RealTensor> & _second_zero;
};

#endif //AUXKERNEL_H
