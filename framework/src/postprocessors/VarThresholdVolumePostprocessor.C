//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VarThresholdVolumePostprocessor.h"

registerMooseObject("MooseApp", VarThresholdVolumePostprocessor);

InputParameters
VarThresholdVolumePostprocessor::validParams()
{
  InputParameters params = VolumePostprocessor::validParams();
  params.addRequiredParam<Real>("threshold",
                                "The value above (or below) which to change the element subdomain");
  params.addParam<MooseEnum>("criterion_type",
                             MooseEnum("BELOW EQUAL ABOVE", "ABOVE"),
                             "Criterion to use for the threshold");
  params.addRequiredCoupledVar("coupled_var",
                               "Coupled variable whose value is used in the criterion");
  params.addClassDescription("Computes the volume of domain that has variable value that is "
                             "below, above, or equal to a certain value.");
  return params;
}

VarThresholdVolumePostprocessor::VarThresholdVolumePostprocessor(const InputParameters & parameters)
  : VolumePostprocessor(parameters),
    _threshold(getParam<Real>("threshold")),
    _criterion_type(getParam<MooseEnum>("criterion_type").getEnum<CriterionType>()),
    _v(coupledValue("coupled_var"))
{
}

void
VarThresholdVolumePostprocessor::execute()
{
  // integrate only if criterion is satisfied
  if (criterionMet() == true)
    ElementIntegralPostprocessor::execute();
}

Real
VarThresholdVolumePostprocessor::computeValue()
{
  Real avg_val = 0;

  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    avg_val += _v[qp] * _JxW[qp] * _coord[qp];
  avg_val /= _current_elem_volume;

  return avg_val;
}

bool
VarThresholdVolumePostprocessor::criterionMet()
{
  auto value = computeValue();
  switch (_criterion_type)
  {
    case CriterionType::Equal:
      return MooseUtils::absoluteFuzzyEqual(value - _threshold, 0);

    case CriterionType::Below:
      return value < _threshold;

    case CriterionType::Above:
      return value > _threshold;

    default:
      mooseError("Unsupported type");
      break;
  }
  return false;
}
