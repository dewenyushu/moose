//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RadialReturnStressUpdate.h"
#include "LinearInterpolation.h"

/**
 * This class uses the stress update material in a radial return isotropic creep
 * model.  This class is one of the basic radial return constitutive models; more complex
 * constitutive models combine creep and plasticity.
 *
 * This class inherits from RadialReturnCreepStressUpdateBase and must be used
 * in conjunction with ComputeMultipleInelasticStress.  This class calculates
 * creep based on stress, temperature, and time effects.  This class also
 * computes the creep strain as a stateful material property.
 */
class ADRateTempDependentStressUpdate : public ADRadialReturnStressUpdate
{
public:
  static InputParameters validParams();

  ADRateTempDependentStressUpdate(const InputParameters & parameters);

  virtual bool substeppingCapabilityEnabled() override;

protected:
  virtual void computeStressInitialize(const ADReal & effective_trial_stress,
                                       const ADRankFourTensor & elasticity_tensor) override;
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & scalar) override;
  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & scalar) override;

  void computePlasticStrainRate(const ADReal & effective_trial_stress,
                         const ADReal & scalar = 0.0);

  void computeMisorientationVariable();
  void computeShearModulusDerivative(const ADRankFourTensor & elasticity_tensor);

  ADReal computeYieldStress();
  ADReal getTemperature();

  virtual void updateInternalStateVariables(const ADReal & effective_trial_stress,
                                            const ADReal & scalar=0.0,
                                            const ADReal & scalar_increment=0.0) override;

  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;
  virtual void computeStressFinalize(const ADRankTwoTensor & plastic_strain_increment) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress, const ADReal & scalar_effective_inelastic_strain) override;

  /// Temperature variable value
  const ADVariableValue * _temperature;

  /// Simulation start time
  const Real _start_time;

  ///
  /// Model parameters (fluid part)
  ///

  /// Melt temperature [K]
  const ADMaterialProperty<Real> & _theta_melt;

  /// Melt viscosity [Pa*s]
  const ADMaterialProperty<Real> & _mu_melt;

  /// Bulk modulus melt [Pa]
  const ADMaterialProperty<Real> & _K_melt;

  ///
  /// Model parameters (solid part)
  ///

  /// Rate independent yield constant [Pa]
  const ADMaterialProperty<Real> & _Y0;

  /// Rate independent yield temperature dependencies [K], [1/K], [K], [-]
  const ADMaterialProperty<Real> & _Y1, & _Y2, & _Y3, & _Y4;

  /// Isotropic hardening shear coefficient [-]
  const ADMaterialProperty<Real> & _Hmu;

  /// Flow rule coefficient constants [1/s], [K]
  const ADMaterialProperty<Real> & _f1, & _f2;

  /// Flow rule exponent constant [-]
  const ADMaterialProperty<Real> & _n1;

  /// Flow rule exponent temperature dependence [K]
  const ADMaterialProperty<Real> & _n2;

  /// Isotropic dynamic recovery constant [Pa]
  const ADMaterialProperty<Real> & _Rd1;

  /// Isotropic dynamic recovery temperature dependence [K]
  const ADMaterialProperty<Real> & _Rd2;

  /// Misorientation variable hardening constant [m/(s Pa)]
  const ADMaterialProperty<Real> & _hxi;

  /// Misorientation variable hardening exponent [-]
  const ADMaterialProperty<Real> & _r;

  /// Components for computing the derivatives
  ADReal _C1, _C2;

  /// Plastic strain rate (flow rule function)
  ADReal _plastic_strain_rate;

  /// Temperature dependent yield stress
  ADReal _yield_stress;

  /// Temperature dependent shear modulus
  ADReal _shear_modulus;

  /// Derivative of shear modulus w.r.t. temperature
  ADReal _shear_modulus_derivative;

  /// LinearInterpolation objects for Young's modulus and Poisson's ratio
  std::unique_ptr<LinearInterpolation> _data_youngs_modulus, _data_poissons_ratio;

  /// Isotropic harderning internal state variable
  ADMaterialProperty<Real> & _hardening_variable;
  const MaterialProperty<Real> & _hardening_variable_old;

  /// Misorientation Variable
  ADMaterialProperty<Real> & _misorientation_variable;
  const MaterialProperty<Real> & _misorientation_variable_old;

  /// Plastic strain material property
  ADMaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
};
