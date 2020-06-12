//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADRateTempDependentStressUpdate.h"

registerMooseObject("TensorMechanicsApp", ADRateTempDependentStressUpdate);

InputParameters
ADRateTempDependentStressUpdate::validParams()
{
  InputParameters params = ADRadialReturnStressUpdate::validParams();
  params.addClassDescription(
      "This class uses the stress update material in a radial return isotropic power law creep "
      "model. This class can be used in conjunction with other creep and plasticity materials "
      "for more complex simulations.");

  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";

  // Linear strain hardening parameters
  params.addCoupledVar("temperature", "Coupled temperature");
  params.addParam<Real>("start_time", 0.0, "Start time (if not zero)");

  // Rate dependent plasticity parameters
  // Default values are from the paper by Michael E. Stender, et. al, 2018.
  params.addParam<Real>("Y0", 5.264e09, "Rate independent yield constant [Pa]");
  params.addParam<Real>("Y1", 2.688e05, "First rate independent yield temperature dependency [K]");
  params.addParam<Real>("Y2", 1.87e-03, "Second rate independent yield temperature dependency [1/K]");
  params.addParam<Real>("Y3", 8.683e02, "Third rate independent yield temperature dependency [K]");
  params.addParam<Real>("Y4", 3.316e01, "Fourth rate independent yield temperature dependency [-]");
  params.addParam<Real>("Hmu", 0.01, "Isotropic hardening shear coefficient [-]");
  params.addParam<Real>("f1", 9.178e-02, "First flow rule coefficient constant [1/s]");
  params.addParam<Real>("f2", 0.0, "Second flow rule coefficient constant [K]");
  params.addParam<Real>("n1", 0.0, "Flow rule exponent constant [-]");
  params.addParam<Real>("n2", 5.699e03, "Flow rule exponent temperature dependence [K]");
  params.addParam<Real>("Rd1", 8.565e02, "Isotropic dynamic recovery constant [Pa]");
  params.addParam<Real>("Rd2", 5.419e03, "Isotropic dynamic recovery temperature dependence [K]");
  params.addParam<Real>("hxi", 1.670e-03, " Misorientation variable hardening constant [m/(s Pa)]");
  return params;
}

ADRateTempDependentStressUpdate::ADRateTempDependentStressUpdate(const InputParameters & parameters)
  : ADRadialReturnStressUpdate(parameters),
    _temperature(isParamValid("temperature") ? &adCoupledValue("temperature") : nullptr),
    _start_time(getParam<Real>("start_time")),
    _Y0(getParam<Real>("Y0")),
    _Y1(getParam<Real>("Y1")),
    _Y2(getParam<Real>("Y2")),
    _Y3(getParam<Real>("Y3")),
    _Y4(getParam<Real>("Y4")),
    _Hmu(getParam<Real>("Hmu")),
    _f1(getParam<Real>("f1")),
    _f2(getParam<Real>("f2")),
    _n1(getParam<Real>("n1")),
    _n2(getParam<Real>("n2")),
    _Rd1(getParam<Real>("Rd1")),
    _Rd2(getParam<Real>("Rd2")),
    _hxi(getParam<Real>("hxi")),
    _plastic_strain(declareADProperty<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "plastic_strain"))
{
  /// initilize _r
  // _r=std::exp(_dG/_G*t);
  // _r=1.0;
}

void
ADRateTempDependentStressUpdate::computeStressInitialize(const ADReal & effective_trial_stress,
                                                     const ADRankFourTensor & elasticity_tensor)
{
  _Y = 0.5*_Y0*(1.0 + std::tanh(_Y2*(_Y3 - (*_temperature)[_qp])))/(_Y4 + std::exp(-_Y1/(*_temperature)[_qp]));

  /// TODO: update temperature dependent shear modulus when temperature changes
  _G = ElasticityTensorTools::getIsotropicShearModulus(elasticity_tensor);
  _dG = -1e-2; // looks like a negative ratio in the plot

  computeFlowRule(effective_trial_stress);

}

ADReal
ADRateTempDependentStressUpdate::computeResidual(const ADReal & /*effective_trial_stress*/,
                                             const ADReal & scalar)
{
  return _phi * _dt - scalar - _C1*_C2*_r*_dG/_G*_dt;
}

ADReal
ADRateTempDependentStressUpdate::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                               const ADReal &  /*scalar*/ )
{
  const ADReal theta = (*_temperature)[_qp];
  const ADReal creep_rate_derivative = _C1*(-3.0*_G/(_r + _Y)) - _C1*_C2*(_Hmu*_G*(1.0+_hxi/_r) -_Rd1*std::exp(-_Rd2/theta)*_r);
  return creep_rate_derivative * _dt - 1.0;
}

void
ADRateTempDependentStressUpdate::computeFlowRule(const ADReal & effective_trial_stress,
                                                      const ADReal & scalar)
{
  const ADReal theta = (*_temperature)[_qp];
  const ADReal stress_delta = effective_trial_stress - 3.0 * _G * scalar;
  const ADReal ratio = stress_delta/(_r + _Y);

  if (ratio>1.0)
  {
    _phi = _f1*std::exp(-_f2*theta)*std::pow(std::sinh(ratio-1.0), _n1+_n2/theta);
    _C1 = (_n1+_n2/theta)*_phi*std::cosh(ratio-1.0)/std::sinh(ratio-1.0);
    _C2 = (3.0*_G*scalar-effective_trial_stress)/(_r + _Y)/(_r + _Y);

    std::cout<<"stress_delta: "<<stress_delta.value()<<" Y:"<<_Y.value()<<" ratio: "<<ratio.value()<<", stress_delta: "<<stress_delta.value()<<", phi: "<<_phi.value()<<"; C1: "<<_C1.value()<<"; C2: "<<_C2.value()<<"; r: "<<_r.value()<<std::endl;
  }
  else
  {
    _phi= 0.0;
    _C1=0.0;
    _C2=0.0;
  }
}

void
ADRateTempDependentStressUpdate::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();

  ADRadialReturnStressUpdate::initQpStatefulProperties();
}

void
ADRateTempDependentStressUpdate::propagateQpStatefulProperties()
{
  _plastic_strain[_qp] = _plastic_strain_old[_qp];

  propagateQpStatefulPropertiesRadialReturn();
}

void
ADRateTempDependentStressUpdate::computeStressFinalize(
    const ADRankTwoTensor & plastic_strain_increment)
{
  _plastic_strain[_qp] = _plastic_strain_old[_qp] + plastic_strain_increment;
}

void
ADRateTempDependentStressUpdate::updateInternalStateVariables(
                                          const ADReal & effective_trial_stress,
                                          const ADReal & scalar)
{
  const ADReal theta = (*_temperature)[_qp];
  ADReal dr = _r*(_dG/_G)+(_Hmu*_G*(1.0+_hxi/_r)-_Rd1*std::exp(-_Rd2/theta)*_r)*scalar;
  _r+=dr;
  computeFlowRule(effective_trial_stress, scalar);
}
