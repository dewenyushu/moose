//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralVectorPostprocessor.h"
#include "libmesh/communicator.h"
#include "Coupleable.h"
#include "MooseVariableDependencyInterface.h"

// Forward Declarations
class MooseMesh;

class AverageValueEveryBlock : public GeneralVectorPostprocessor,
                               public Coupleable,
                               public MooseVariableDependencyInterface
{
public:
  static InputParameters validParams();

  AverageValueEveryBlock(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

protected:
  // Compute the area as well as integral of the coupled variables on the current element
  void computeIntegral();

  /// Reference to the mesh
  MooseMesh & _mesh;
  Assembly & _assembly;
  const MooseArray<Point> & _q_point;
  const QBase * const & _qrule;
  const MooseArray<Real> & _JxW;
  const MooseArray<Real> & _coord;

  /// Variables to output
  std::vector<VariableName> _variables;

  std::vector<const VariableValue *> _variable_vals;

  /// Number of columns
  int _num_cols;

  /// Number of rows
  int _num_rows;

  /// Area for each subdomain
  std::vector<Real> _subdomain_areas;

  /// Vector of outputs, where each entry is the vector of average values for single variable in each block
  std::vector<VectorPostprocessorValue *> _output_vector;

  /// Map between the row index and subdomain ID
  std::unordered_map<SubdomainID, int> _block_id_map;
};
