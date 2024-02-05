//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#ifdef NEML2_ENABLED
#include "neml2/tensors/LabeledVector.h"
#include "neml2/tensors/LabeledMatrix.h"
#include "neml2/models/Model.h"
#endif

#include "NEML2SolidMechanicsInterface.h"

#include "BatchMaterial.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "SymmetricRankTwoTensor.h"
#include "SymmetricRankFourTensor.h"

#include "BatchScalarProperty.h"

typedef BatchMaterial<BatchMaterialUtils::TupleStd,
                      // Outputs: stress, internal variables, dstress/dstrain, dstress/dparam
                      std::tuple<RankTwoTensor, RankFourTensor, RankTwoTensor>,
                      // Inputs:
                      //   strain
                      //   temperature
                      BatchMaterialUtils::GatherMatProp<RankTwoTensor>,
                      BatchMaterialUtils::GatherVariable>
    CauchyStressFromNEML2UOParent;

/**
 * This user object gathers input variables required for an objective stress integration from all
 * quadrature points. The batched input vector is sent through a NEML2 material model to perform the
 * constitutive update.
 */
class CauchyStressFromNEML2UO : public NEML2SolidMechanicsInterface<CauchyStressFromNEML2UOParent>
{
public:
  static InputParameters validParams();

  CauchyStressFromNEML2UO(const InputParameters & params);

#ifndef NEML2_ENABLED
  virtual void preCompute() {}
  virtual void batchCompute() override {}
  virtual void postCompute() {}
#else

  virtual void batchCompute() override;
  virtual void preCompute();
  virtual void timestepSetup() override;
  virtual void postCompute();

protected:
  /// Advance state and forces in time
  virtual void advanceStep();

  /// Update the forces driving the material model update
  virtual void updateForces();

  /// Apply the predictor to set current trial state
  virtual void applyPredictor();

  /// Perform the material update
  virtual void solve();

  /// The input vector of the material model
  neml2::LabeledVector _in;

  /// The output vector of the material model
  neml2::LabeledVector _out;

  /// The derivative of the output vector w.r.t. the input vector
  neml2::LabeledMatrix _dout_din;

  /// List of model parameters for which we wish to compute derivatives for
  const std::vector<std::string> & _parameter_derivatives;

  /// Flag to check whether derivative w.r.t. model parameters are requested
  const bool _require_parameter_derivatives;

  /// List of model parameters for which we wish set from MOOSE batch materials
  const std::vector<std::string> & _reset_parameter_names;

  /// Flag to check whether we wish to set model parameters from MOOSE batch materials
  const bool _reset_parameters;

  const MaterialProperty<Real> * _material_param;

  const BatchScalarProperty * _material_param_uo;

#endif // NEML2_ENABLED
};
