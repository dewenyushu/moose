//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VolumePostprocessor.h"

/**
 * This postprocessor computes the volume of a specified block.
 */
class VarThresholdVolumePostprocessor : public VolumePostprocessor
{
public:
  static InputParameters validParams();

  VarThresholdVolumePostprocessor(const InputParameters & parameters);

  virtual void execute() override;

private:
  Real computeValue();

  bool criterionMet();

  /// Threshold to modify the element subdomain ID
  const Real _threshold;

  /// Criterion type
  const enum class CriterionType { Below, Equal, Above } _criterion_type;

  /// The coupled variable used in the criterion
  const VariableValue & _v;
};
