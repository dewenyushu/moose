//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VarCenterGaussianHeatSource.h"

registerMooseObject("CombinedApp", VarCenterGaussianHeatSource);

InputParameters
VarCenterGaussianHeatSource::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<Real>("power", "power");
  params.addParam<Real>("efficiency", 1, "process efficiency");
  params.addParam<bool>("use_input_r",
                        true,
                        "option to use user input effective radii or from experimentally fitted "
                        "formulations. Default is to use user input data.");
  params.addParam<bool>("auto_heat_source_height",
                        false,
                        "option to detect the current material height and use it as the "
                        "z-coordinate of the heat source center.");
  params.addParam<std::vector<Real>>(
      "r",
      {},
      "effective radii along three directions. If only one parameter is provided, then we assume "
      "the effective radius to be equal along three directions.");

  params.addParam<Real>(
      "feed_rate",
      0.000124,
      "powder material feed rate. This value is used only when use_input_r = false.");

  params.addParam<Real>(
      "factor",
      1.0,
      "scaling factor that is multiplied to the heat source to adjust the intensity");

  params.addParam<Real>("std_factor", 1.0, "factor to the std to sample r");

  params.addRequiredCoupledVar(
      "position_x", "Coupled variable representing the x-coordinate of the heat source center.");

  params.addRequiredCoupledVar(
      "position_y", "Coupled variable representing the y-coordinate of the heat source center.");

  params.addRequiredCoupledVar(
      "position_z", "Coupled variable representing the z-coordinate of the heat source center.");

  params.addCoupledVar(
      "speed_x", "Coupled variable representing the x-direction speed of the heat source center.");

  params.addCoupledVar(
      "speed_y", "Coupled variable representing the y-direction speed of the heat source center.");

  params.addCoupledVar(
      "speed_z", "Coupled variable representing the z-direction speed of the heat source center.");

  params.addParam<std::string>("material_height_uo",
                               "name of the user object that tracks the material height.");
  params.declareControllable("power efficiency r factor");

  params.addClassDescription(
      "Gaussian volumetric source heat source using center represented by (aux) variables.");

  return params;
}

VarCenterGaussianHeatSource::VarCenterGaussianHeatSource(const InputParameters & parameters)
  : Material(parameters),
    _use_input_r(getParam<bool>("use_input_r")),
    _auto_heat_source_height(getParam<bool>("auto_heat_source_height")),
    _P(getParam<Real>("power")),
    _eta(getParam<Real>("efficiency")),
    _feed_rate(getParam<Real>("feed_rate")),
    _r(getParam<std::vector<Real>>("r")),
    _f(getParam<Real>("factor")),
    _std_factor(getParam<Real>("std_factor")),
    _volumetric_heat(declareADProperty<Real>("volumetric_heat")),
    _material_height_uo_name(
        isParamValid("material_height_uo") ? getParam<std::string>("material_height_uo") : ""),
    _material_height_uo(
        isParamValid("material_height_uo")
            ? &_fe_problem.getUserObject<GeometryHeightTracker>(_material_height_uo_name)
            : nullptr),
    _px(coupledValue("position_x")),
    _py(coupledValue("position_y")),
    _pz(coupledValue("position_z")),
    _vx(isParamValid("speed_x") ? &coupledValue("speed_x") : nullptr),
    _vy(isParamValid("speed_y") ? &coupledValue("speed_y") : nullptr),
    _vz(isParamValid("speed_z") ? &coupledValue("speed_z") : nullptr)
{
  if (_use_input_r)
  {
    if (_r.size() != 1 && _r.size() != 3)
      paramError("r", "The effective radii should have 1 or 3 components.");
    // make sure we have 3 equal components if only one parameter is provided
    if (_r.size() == 1)
    {
      _r.push_back(_r[0]);
      _r.push_back(_r[0]);
    }
  }
  else
  {
    _r.resize(3);
    // make sure we have the speed vars in this case
    if (!isParamValid("speed_x"))
      paramError("speed_x", "The speed_x parameter needs to be valid when use_input_r = false");
    if (!isParamValid("speed_y"))
      paramError("speed_y", "The speed_y parameter needs to be valid when use_input_r = false");
    if (!isParamValid("speed_z"))
      paramError("speed_z", "The speed_z parameter needs to be valid when use_input_r = false");
  }

  if (_auto_heat_source_height && !isParamValid("material_height_uo"))
    paramError("material_height_uo",
               "A valid material height user object is needed when auto_heat_source_height=true.");

  // set to a small number to start
  _r_time_prev = -1.0e8;
}

void
VarCenterGaussianHeatSource::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  _volumetric_heat[_qp] = computeHeatSource(x, y, z);
}

Real
VarCenterGaussianHeatSource::computeHeatSource(const Real x, const Real y, const Real z)
{
  // center of the heat source
  Real x_t = _px[_qp];
  Real y_t = _py[_qp];
  Real z_t = _pz[_qp];

  // overwrite the z_t if we automatically track the height of material
  if (_auto_heat_source_height)
    z_t = _material_height_uo->getLayerHeight();

  Real dist_x = -2.0 * std::pow(x - x_t, 2.0) / _r[0] / _r[0];
  Real dist_y = -2.0 * std::pow(y - y_t, 2.0) / _r[1] / _r[1];
  Real dist_z = -2.0 * std::pow(z - z_t, 2.0) / _r[2] / _r[2];

  // effective radii under current processing parameters
  if (!_use_input_r && _t > _r_time_prev)
  {
    // compute scanning speed at this time
    _scan_speed = std::sqrt((*_vx)[_qp] * (*_vx)[_qp] + (*_vy)[_qp] * (*_vy)[_qp] +
                            (*_vz)[_qp] * (*_vz)[_qp]);
    // compute the effective radii
    computeEffectiveRadii();

    // update time
    _r_time_prev = _t;
  }

  // Gaussian point heat source
  return 2.0 * _P * _eta * _f / libMesh::pi / _r[0] / _r[1] / _r[2] *
         std::exp(dist_x + dist_y + dist_z);
}

void
VarCenterGaussianHeatSource::computeEffectiveRadii()
{
  // get scaled laser power, scanning speed, and powder feed rate
  Real lp = _P / 1.0e-3 / 400.0;             // input in unit 1e-3 W
  Real ss = _scan_speed / 4.23333e-4 / 40.0; // input in mm/ms (1 ipm = 4.23333e-4 mm/ms)
  Real pf = _feed_rate / 0.031e-3 / 15.0;    // input in g/ms (1rpm = 0.031e-3 g/ms)

  // list of the variable values
  std::vector<Real> vals = {lp, ss, pf, lp * lp, ss * ss, pf * pf, lp * ss, lp * pf, ss * pf, 1.0};

  Real mean_rxy = 0.0, mean_rz = 0.0;
  Real std_r = 0.0;
  _r[2] = 0.0;
  for (unsigned int i = 0; i < vals.size(); ++i)
  {
    // compute the mean and standard deviation of the melt pool dimension (x and y dimensions)
    mean_rxy += _diameter_param[i].first * vals[i];
    std_r += _diameter_param[i].second * vals[i];
    // compute the material bead height (z dimension)
    mean_rz += _height_param[i] * vals[i];
  }

  // sample the r[0] and r[2] value
  // define a random number generator
  std::random_device rd{};
  std::mt19937 generator{rd()};
  std::normal_distribution<double> dist_xy(mean_rxy, std_r);
  std::normal_distribution<double> dist_z(mean_rz, std_r);
  _r[0] = dist_xy(generator);
  _r[2] = dist_z(generator);
  // make sure that the sampled values are within +-sigma range
  if (_r[0] > mean_rxy + _std_factor * std::abs(std_r))
    _r[0] = mean_rxy + _std_factor * std::abs(std_r);
  else if (_r[0] < mean_rxy - _std_factor * std::abs(std_r) || _r[0] <= 0.0)
    _r[0] = std::abs(mean_rxy - _std_factor * std::abs(std_r)); // make sure r is not negative
  _r[0] *= 0.5;                                                 // get the radius
  _r[1] = _r[0];

  if (_r[2] > mean_rz + _std_factor * std::abs(std_r))
    _r[2] = mean_rz + _std_factor * std::abs(std_r);
  else if (_r[2] < mean_rz - _std_factor * std::abs(std_r) || _r[0] <= 0.0)
    _r[2] = std::abs(mean_rz - _std_factor * std::abs(std_r)); // make sure r is not negative
}
