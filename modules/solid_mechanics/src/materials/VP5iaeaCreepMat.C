//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VP5iaeaCreepMat.h"

registerMooseObject("SolidMechanicsApp", VP5iaeaCreepMat);
registerMooseObject("SolidMechanicsApp", ADVP5iaeaCreepMat);

template <bool is_ad>
InputParameters
VP5iaeaCreepMatTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription(
      "This class uses the stress update material in a radial return isotropic power law creep "
      "model. This class can be used in conjunction with other creep and plasticity materials "
      "for more complex simulations.");

  // Creep strain parameters
  params.addRequiredCoupledVar("dose", "Coupled damage dose");
  params.addRequiredParam<Real>("youngs_modulus", "Youngs modulus");
  return params;
}

template <bool is_ad>
VP5iaeaCreepMatTempl<is_ad>::VP5iaeaCreepMatTempl(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _dose(this->coupledValue("dose")),
    _dose_old(this->coupledValueOld("dose")),
    _youngs_modulus(this->template getParam<Real>("youngs_modulus"))
{
}

template <bool is_ad>
template <typename ScalarType>
ScalarType
VP5iaeaCreepMatTempl<is_ad>::computeResidualInternal(
    const GenericReal<is_ad> & effective_trial_stress, const ScalarType & scalar)
{
  const ScalarType stress_delta =
      effective_trial_stress - _three_shear_modulus * scalar; // <-- this is not used
  const ScalarType creep_increment =
      stress_delta / _youngs_modulus * (_dose[_qp] - _dose_old[_qp]) / 4.0;
  return creep_increment - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
VP5iaeaCreepMatTempl<is_ad>::computeDerivative(const GenericReal<is_ad> & /* effective_trial_stress */,
                                               const GenericReal<is_ad> & /* scalar */)
{
  const GenericReal<is_ad> creep_increment_derivative =
      -_three_shear_modulus / _youngs_modulus * (_dose[_qp] - _dose_old[_qp]) / 4.0;
  return creep_increment_derivative - 1.0;
}

template <bool is_ad>
void
VP5iaeaCreepMatTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
VP5iaeaCreepMatTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
VP5iaeaCreepMatTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->_use_substepping != RadialReturnStressUpdateTempl<is_ad>::SubsteppingType::NONE;
}

template class VP5iaeaCreepMatTempl<false>;
template class VP5iaeaCreepMatTempl<true>;
template Real VP5iaeaCreepMatTempl<false>::computeResidualInternal<Real>(const Real &,
                                                                         const Real &);
template ADReal VP5iaeaCreepMatTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                            const ADReal &);
template ChainedReal
VP5iaeaCreepMatTempl<false>::computeResidualInternal<ChainedReal>(const Real &,
                                                                  const ChainedReal &);
template ChainedADReal
VP5iaeaCreepMatTempl<true>::computeResidualInternal<ChainedADReal>(const ADReal &,
                                                                   const ChainedADReal &);
