//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "MoosePreconditioner.h"

// libMesh includes
#include "libmesh/preconditioner.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/mesh_base.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <vector>

// Forward declarations
class NonlinearSystemBase;
class ParallelDualMortarPreconditioner;

template <>
InputParameters validParams<ParallelDualMortarPreconditioner>();

/**
 * Interface for condensing out LMs for the dual mortar approach.
 */
class ParallelDualMortarPreconditioner : public MoosePreconditioner, public Preconditioner<Number>
{
public:
  static InputParameters validParams();

  ParallelDualMortarPreconditioner(const InputParameters & params);
  virtual ~ParallelDualMortarPreconditioner();

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init();

  /**
   * This is called every time the "operator might have changed".
   *
   * This is essentially where you need to fill in your preconditioning matrix.
   */
  virtual void setup();

  /**
   * Computes the preconditioned vector "x" based on input "y".
   * Usually by solving Px=y to get the action of P^-1 y.
   */
  virtual void apply(const NumericVector<Number> & y, NumericVector<Number> & x);

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear();

protected:
  /**
   * Get dofs for the variable to be condensed out
   */
  void getDofToCondense();

  /**
   * Get dofs for the primary variable on the contact interface
   */
  void getDofContact();

  /**
   * Get row and col dofs for the condensed system
   */
  void getDofColRow();

  /**
   * Reconstruct the equation system
   */
  void condenseSystem();

  /**
   * Get condensed x and y
   */
  void getCondensedXY(const NumericVector<Number> & y, NumericVector<Number> & x);

  /**
   * Compute Lagrange multipliers using updated solution vector
   */
  void computeLM();

  /**
   * Assemble the full solution vector
   */
  void getFullSolution(const NumericVector<Number> & y, NumericVector<Number> & x);

  /**
   * Print active node info for debugging purposes
   */
  void print_node_info();

  /// The nonlinear system this PC is associated with (convenience reference)
  NonlinearSystemBase & _nl;
  /// Mesh object for easy reference
  MooseMesh * _mesh;
  /// DofMap for easy reference
  DofMap * _dofmap;
  /// Number of variables
  unsigned int _n_vars;
  /// Name and ID of the variable that is to be condensed out (usually the Lagrange multiplier variable)
  const std::vector<std::string> _var_names;
  std::vector<unsigned int> _var_ids;
  // Name and ID of the corresponding coupled variable
  const std::vector<std::string> _cp_var_names;
  std::vector<unsigned int> _cp_var_ids;

  /// Submatrices (_1 -> primary subdomain; _2 -> secondary subdomain; _i -> interior; _c -> contact interface)
  std::unique_ptr<PetscMatrix<Number>> _D, _M, _MDinv, _u2c_rows;

  /// Condensed Jacobian
  std::unique_ptr<PetscMatrix<Number>> _J_condensed;

  /// Condensed x and y
  std::unique_ptr<NumericVector<Number>> _x_hat, _y_hat, _r2c, _lambda;

  /// Contact
  const BoundaryID _primary_boundary, _secondary_boundary;

  /// Subdomain
  const SubdomainID _primary_subdomain, _secondary_subdomain;

  /// Whether DOFs info has been saved
  mutable bool _save_dofs;

  /// Which preconditioner to use for the solve.
  PreconditionerType _pre_type;

  /// Holds one Preconditioner object per small system to solve.
  std::unique_ptr<Preconditioner<Number>> _preconditioner;

  /// indices
  std::vector<numeric_index_type> _glm, _lm, _gu1c, _u1c, _gu2c, _u2c;

  std::vector<numeric_index_type> _grows, _rows, _gcols, _cols;

  /// map of glm to gu2c
  std::map<numeric_index_type, numeric_index_type> _map_glm_gu2c, _map_gu2c_glm;

  // gu2c index as keys, the corresponding row index in _D as the value
  std::map<numeric_index_type, numeric_index_type> _map_gu2c_to_order;

  /// Timers
  PerfID _init_timer;
  PerfID _apply_timer;
};
