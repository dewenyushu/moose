//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeCrystalPlasticityThermalEigenstrain.h"

registerMooseObject("TensorMechanicsApp", ComputeCrystalPlasticityThermalEigenstrain);

InputParameters
ComputeCrystalPlasticityThermalEigenstrain::validParams()
{
  InputParameters params = ComputeCrystalPlasticityEigenstrainBase::validParams();

  params.addCoupledVar("temperature", "Coupled temperature variable");
  params.addParam<std::vector<Real>>(
      "thermal_expansion_coefficients",
      "Vector of values defining the constant second order thermal expansion coefficients");

  return params;
}

ComputeCrystalPlasticityThermalEigenstrain::ComputeCrystalPlasticityThermalEigenstrain(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ComputeCrystalPlasticityEigenstrainBase>(parameters),
    _temperature(coupledValue("temperature")),
    _temperature_old(coupledValueOld("temperature")),
    _ddeformation_gradient_dT(declarePropertyDerivative<RankTwoTensor>(
        _deformation_gradient_name, getVar("temperature", 0)->name()))
{
  std::vector<Real> therm_exp_coeff = getParam<std::vector<Real>>("thermal_expansion_coefficients");
  if (therm_exp_coeff.size() != _mesh.dimension())
    paramError("thermal_expansion_coefficients",
               "thermal expansion coefficient size must match the mesh dimension");
  _thermal_expansion_coefficients.fillFromInputVector(therm_exp_coeff);
}

void
ComputeCrystalPlasticityThermalEigenstrain::computeQpDeformationGradient()
{
  // compute the deformation gradient due to thermal expansion
  Real dtheta = (_temperature[_qp] - _temperature_old[_qp]) * _substep_dt / _dt;
  RankTwoTensor residual_equivalent_thermal_expansion_increment =
      RankTwoTensor::Identity() - dtheta * _thermal_expansion_coefficients;
  RankTwoTensor _inverse_deformation_gradient =
      _deformation_gradient_old[_qp].inverse() * residual_equivalent_thermal_expansion_increment;
  _deformation_gradient[_qp] = _inverse_deformation_gradient.inverse();

  // compute the derivative of deformation gradient w.r.t temperature
  _ddeformation_gradient_dT[_qp] = _thermal_expansion_coefficients * _deformation_gradient[_qp];
}
