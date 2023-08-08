//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "UraniumDioxideCrystalPlasticityUpdate.h"
#include "DelimitedFileReader.h"

registerMooseObject("TensorMechanicsApp", UraniumDioxideCrystalPlasticityUpdate);

InputParameters
UraniumDioxideCrystalPlasticityUpdate::validParams()
{
  InputParameters params = CrystalPlasticityStressUpdateBase::validParams();
  params.addClassDescription("Crystal plasticity model for UO2.");
  params.addRequiredParam<FileName>(
      "parameter_filename", "file name of the text file providing the slip resistance information");
  params.addParam<Real>("dG0",
                        4,
                        "effective thermal activation energy without any external stress for "
                        "dislocation glide, 3 - 7 [eV]");
  params.addParam<Real>("k", 8.617e-5, "Boltzmann constant [eV/K]");
  params.addParam<Real>("T", 300, "absolute temperature [K]");
  params.addParam<Real>(
      "p",
      0.7,
      "[0-1], 1st exponent parameter related to the shape of the obstacles resistance profile");
  params.addParam<Real>(
      "q",
      1.4,
      "[1-2], 2nd exponent parameter related to the shape of the obstacles resistance profile");
  params.addParam<unsigned int>("nu",
                                1.2e10,
                                "attempt (or attack) frequency, which is the number of counts that "
                                "a dislocation attempts to overcome an obstacle [s^-1]");
  params.addParam<Real>("lambda_eff", 1.0, "mean distance between obstacles");
  params.addParam<Real>("tau0", 20, "lattice friction or Peierls stress, 20-50 [MPa]");
  params.addParam<Real>("mu", 1e5, "shear modulus [MPa]");
  params.addParam<Real>("b", 3.87e-7, "magnitude of the burger's vector [mm]");
  params.addParam<Real>(
      "k1", 2e6, "constant in the dislocation density evolutino equation [mm^-1]");
  params.addParam<Real>("chi", 0.5, "interaction parameter, 0.1-0.9");
  params.addParam<Real>("gamma_dot0", 1e7, "the reference shear rate [s^-1]");
  params.addParam<Real>("rho0", 1e-15, "initial forest dislocation density value");

  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");

  return params;
}

UraniumDioxideCrystalPlasticityUpdate::UraniumDioxideCrystalPlasticityUpdate(
    const InputParameters & parameters)
  : CrystalPlasticityStressUpdateBase(parameters),
    // filename containing all the pertaining parameters
    _parameter_filename(getParam<FileName>("parameter_filename")),
    // Parameters in the constitutive equations
    _dG0(getParam<Real>("dG0")),
    _k(getParam<Real>("k")),
    _T(getParam<Real>("T")),
    _p(getParam<Real>("p")),
    _q(getParam<Real>("q")),
    _nu(getParam<unsigned int>("nu")),
    _lambda_eff(getParam<Real>("lambda_eff")),
    _tau0(getParam<Real>("tau0")),
    _mu(getParam<Real>("mu")),
    _b(getParam<Real>("b")),
    _k1(getParam<Real>("k1")),
    _chi(getParam<Real>("chi")),
    _gamma_dot0(getParam<Real>("gamma_dot0")),
    _rho0(getParam<Real>("rho0")),
    // state variables
    _forest_dislocation_density(
        declareProperty<std::vector<Real>>(_base_name + "forest_dislocation_density")),
    _forest_dislocation_density_old(
        getMaterialPropertyOld<std::vector<Real>>(_base_name + "forest_dislocation_density")),

    _slip_resistance_increment(_number_slip_systems, 0.0),
    // resize local caching vectors used for substepping
    _previous_substep_slip_resistance(_number_slip_systems, 0.0),
    _slip_resistance_before_update(_number_slip_systems, 0.0)
{
  // get constitutive parameters from file
  readDislocationSlipResistanceFileParams();

  // check range of certain parameters
  if (_p < 0.0 || _p > 1.0)
    paramError("p", "the 1st exponent parameter should range between 0 and 1.");
  if (_q < 1.0 || _q > 2.0)
    paramError("q", "the 2nd exponent parameter should range between 1 and 2.");
}

void
UraniumDioxideCrystalPlasticityUpdate::readDislocationSlipResistanceFileParams()
{
  _g.resize(_number_slip_systems);
  _D.resize(_number_slip_systems);
  _K.resize(_number_slip_systems);

  // Read the data from CSV file
  MooseUtils::DelimitedFileReader reader(_parameter_filename);
  // data in rows
  reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  reader.read();

  // number of columns for properties including activation energy,
  // drag stress, and the hardening matrix
  unsigned int num_cols = 2 + _number_slip_systems;

  const std::vector<std::vector<double>> & data = reader.getData();
  if (data.size() != _number_slip_systems)
    paramError("parameter_filename", "number of rows incorrect in parameter info file");
  if (data[0].size() != num_cols)
    paramError("parameter_filename", "number of columns incorrect in parameter info file");

  // start grabbing the data
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
    _g[i] = data[i][0];
    _D[i] = data[i][1];
    // K should be a square matrix considering latent/self hardening
    _K[i] = std::vector<Real>(data[i].begin() + 2, data[i].end());
  }
}

void
UraniumDioxideCrystalPlasticityUpdate::initQpStatefulProperties()
{
  CrystalPlasticityStressUpdateBase::initQpStatefulProperties();

  _forest_dislocation_density[_qp].resize(_number_slip_systems);
  _forest_dislocation_density_rate.resize(_number_slip_systems);

  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
    _forest_dislocation_density[_qp][i] = _rho0;
    _slip_resistance[_qp][i] = _tau0;
    _slip_increment[_qp][i] = _gamma_dot0;
  }
}

void
UraniumDioxideCrystalPlasticityUpdate::setInitialConstitutiveVariableValues()
{
  // Would also set old dislocation densities here if included in this model
  _slip_resistance[_qp] = _slip_resistance_old[_qp];
  _previous_substep_slip_resistance = _slip_resistance_old[_qp];

  _forest_dislocation_density[_qp] = _forest_dislocation_density_old[_qp];
  _previous_substep_forest_dislocation_density = _forest_dislocation_density_old[_qp];
}

void
UraniumDioxideCrystalPlasticityUpdate::setSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _slip_resistance[_qp] = _previous_substep_slip_resistance;
  _forest_dislocation_density[_qp] = _previous_substep_forest_dislocation_density;
}

bool
UraniumDioxideCrystalPlasticityUpdate::calculateSlipRate()
{
  // follow the flow rule in equation 2
  for (const auto i : make_range(_number_slip_systems))
  {
    auto stress_effect = 1.0 - std::pow(std::abs(_tau[_qp][i]) / _slip_resistance[_qp][i], _p);
    auto temp_stress_effects = -_dG0 / _k / _T * std::pow(stress_effect, _q);
    _slip_increment[_qp][i] = _forest_dislocation_density[_qp][i] * _b * _lambda_eff * _nu *
                              std::exp(temp_stress_effects);

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {
      if (_print_convergence_message)
        mooseWarning("Maximum allowable slip increment exceeded ",
                     std::abs(_slip_increment[_qp][i]) * _substep_dt);

      return false;
    }
  }
  return true;
}

void
UraniumDioxideCrystalPlasticityUpdate::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  for (const auto i : make_range(_number_slip_systems))
  {
    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0))
      dslip_dtau[i] = 0.0; // derivative is actually not defined if shear stress equals 0
    else
    {
      auto stress_effect = 1.0 - std::pow(std::abs(_tau[_qp][i]) / _slip_resistance[_qp][i], _p);
      auto stress_effect_deriv =
          _p * std::pow(std::abs(_tau[_qp][i]) / _slip_resistance[_qp][i], _p - 1);
      // sign change due to derivative of the absolute value of the shear stress
      if (_tau[_qp][i] < 0.0)
        stress_effect_deriv *= -1.0;
      auto temp_stress_deriv =
          -_dG0 / _k / _T * _q * std::pow(stress_effect, _q - 1.0) * stress_effect_deriv;
      dslip_dtau[i] = _slip_increment[_qp][i] * temp_stress_deriv;
    }
  }
}

bool
UraniumDioxideCrystalPlasticityUpdate::areConstitutiveStateVariablesConverged()
{
  return isConstitutiveStateVariableConverged(_slip_resistance[_qp],
                                              _slip_resistance_before_update,
                                              _previous_substep_slip_resistance,
                                              _resistance_tol);
}

void
UraniumDioxideCrystalPlasticityUpdate::updateSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _previous_substep_slip_resistance = _slip_resistance[_qp];

  _previous_substep_forest_dislocation_density = _forest_dislocation_density[_qp];
}

void
UraniumDioxideCrystalPlasticityUpdate::cacheStateVariablesBeforeUpdate()
{
  _slip_resistance_before_update = _slip_resistance[_qp];

  _forest_dislocation_density_before_update = _forest_dislocation_density[_qp];
}

void
UraniumDioxideCrystalPlasticityUpdate::calculateStateVariableEvolutionRateComponent()
{
  /*
   * Calculate the forest dislocation density rate.
   */
  for (const auto i : make_range(_number_slip_systems))
  {
    // calculate k2 for i-th slip system (follow equation 10)
    auto k2 = _k1 * _chi * _b / _g[i] *
              (1.0 - _k * _T / (_D[i] * std::pow(_b, 3)) *
                         std::log(_slip_increment[_qp][i] / _gamma_dot0));
    // calculate forest dislocation density for i-th slip system (follow equation 9)
    _forest_dislocation_density_rate[i] =
        _slip_increment[_qp][i] * (_k1 * std::sqrt(_forest_dislocation_density[_qp][i]) -
                                   k2 * _forest_dislocation_density[_qp][i]);
  }
}

bool
UraniumDioxideCrystalPlasticityUpdate::updateStateVariables()
{
  /*
   * Update the forest dislocation density using the rate and previous substep value.
   * Negative dislocation density indicates unconverged solution.
   */
  for (const auto i : make_range(_number_slip_systems))
  {
    auto dislocation_density_increment = _forest_dislocation_density_rate[i] * _substep_dt;

    if (_previous_substep_forest_dislocation_density[i] < 0.0 &&
        _forest_dislocation_density_rate[i] < 0.0)
      _forest_dislocation_density[_qp][i] = _previous_substep_forest_dislocation_density[i];
    else
      _forest_dislocation_density[_qp][i] =
          _previous_substep_forest_dislocation_density[i] + dislocation_density_increment;

    if (_forest_dislocation_density[_qp][i] < 0.0)
      return false;
  }
  return true;
}

void
UraniumDioxideCrystalPlasticityUpdate::calculateSlipResistance()
{
  /*
   * Update the slip resistance, considering forest dislocations only.
   */
  for (const auto i : make_range(_number_slip_systems))
  {
    Real sum = 0.0;
    for (const auto j : make_range(_number_slip_systems))
      sum += _K[i][j] * _forest_dislocation_density[_qp][j];
    // Note: forest dislocation resistance can be a material property as well, if needed in the
    // phase-field model in the future
    auto forest_dislocation_resistance = _mu * _b * std::sqrt(sum);

    _slip_resistance[_qp][i] = _tau0 + std::abs(forest_dislocation_resistance);
  }
}
