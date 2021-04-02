//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeCrystalPlasticityEigenstrainBase.h"
#include "RankTwoTensor.h"

InputParameters
ComputeCrystalPlasticityEigenstrainBase::validParams()
{
  InputParameters params = Material::validParams();

  // The return stress increment classes are intended to be iterative materials, so must set compute
  // = false for all inheriting classes
  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addRequiredParam<std::string>("eigenstrain_name",
                                       "Material property name for the eigenstrain tensor computed "
                                       "by this model.");
  params.addRequiredParam<std::string>(
      "deformation_gradient_name",
      "Material property name for the deformation gradient tensor computed "
      "by this model.");
  return params;
}

ComputeCrystalPlasticityEigenstrainBase::ComputeCrystalPlasticityEigenstrainBase(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _eigenstrain_name(_base_name + getParam<std::string>("eigenstrain_name")),
    _deformation_gradient_name(_base_name + getParam<std::string>("deformation_gradient_name")),
    _eigenstrain(declareProperty<RankTwoTensor>(_eigenstrain_name)),
    _deformation_gradient(declareProperty<RankTwoTensor>(_deformation_gradient_name)),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>(_deformation_gradient_name)),
    _step_zero(declareRestartableData<bool>("step_zero", true))
{
}

void
ComputeCrystalPlasticityEigenstrainBase::initQpStatefulProperties()
{
  // This property can be promoted to be stateful by other models that use it,
  // so it needs to be initalized.
  _eigenstrain[_qp].zero();

  // Initialize deformation gradient to be identity
  _deformation_gradient[_qp].zero();
  _deformation_gradient[_qp].addIa(1.0);
}

void
ComputeCrystalPlasticityEigenstrainBase::computeQpProperties()
{
  if (_t_step >= 1)
    _step_zero = false;

  // Skip the eigenstrain calculation in step zero because no solution is computed during
  // the zeroth step, hence computing the eigenstrain in the zeroth step would result in
  // an incorrect calculation of mechanical_strain, which is stateful.
  if (!_step_zero)
  {
    computeQpDeformationGradient();
    computeQpEigenstrain();
  }
}

void
ComputeCrystalPlasticityEigenstrainBase::computeQpEigenstrain()
{
  _eigenstrain[_qp] = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] -
                      RankTwoTensor::Identity();
  _eigenstrain[_qp] = _eigenstrain[_qp] * 0.5;
}

void
ComputeCrystalPlasticityEigenstrainBase::setQp(const unsigned int & qp)
{
  _qp = qp;
}

void
ComputeCrystalPlasticityEigenstrainBase::setSubstepDt(const Real & substep_dt)
{
  _substep_dt = substep_dt;
}

const RankTwoTensor
ComputeCrystalPlasticityEigenstrainBase::getDeformationGradient()
{
  return _deformation_gradient[_qp];
}

const RankTwoTensor
ComputeCrystalPlasticityEigenstrainBase::getDeformationGradientInverse()
{
  return _deformation_gradient[_qp].inverse();
}
