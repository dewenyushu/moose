//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FunctionPathEllipsoidHeatSource.h"

registerMooseObject("HeatConductionApp", FunctionPathEllipsoidHeatSource);

template <>
InputParameters
validParams<FunctionPathEllipsoidHeatSource>()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<Real>("power", "power");
  params.addParam<Real>("efficiency", 1, "process efficiency");
  params.addRequiredParam<Real>("r", "effective radius");
  params.addParam<Real>(
      "factor", 1, "scaling factor that is multiplied to the heat source to adjust the intensity");
  params.addParam<FunctionName>(
      "function_x", "0", "The x component of the center of the heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_y", "0", "The y component of the center of the heating spot as a function of time");
  params.addParam<FunctionName>(
      "function_z", "0", "The z component of the center of the heating spot as a function of time");

  params.addParam<MooseEnum>(
      "heat_source_type", MooseEnum("point line mixed", "point"), "Type of the heat source");

  params.addParam<Real>(
      "threshold_length", 1.0, "Threshold size when we change the way of computing heat source");

  params.addParam<Real>(
      "number_time_integration",
      10,
      "Number of points to do time integration for averaged heat source calculation");

  params.addClassDescription("Double ellipsoid volumetric source heat with function path.");

  return params;
}

FunctionPathEllipsoidHeatSource::FunctionPathEllipsoidHeatSource(const InputParameters & parameters)
  : Material(parameters),
    _P(getParam<Real>("power")),
    _eta(getParam<Real>("efficiency")),
    _r(getParam<Real>("r")),
    _f(getParam<Real>("factor")),
    _function_x(getFunction("function_x")),
    _function_y(getFunction("function_y")),
    _function_z(getFunction("function_z")),
    _heat_source_type(getParam<MooseEnum>("heat_source_type").getEnum<HeatSourceType>()),
    _threshold_length(getParam<Real>("threshold_length")),
    _number_time_integration(getParam<Real>("number_time_integration")),
    _volumetric_heat(declareADProperty<Real>("volumetric_heat"))
{
}

void
FunctionPathEllipsoidHeatSource::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  switch (_heat_source_type)
  {
    case HeatSourceType::POINT:
      _volumetric_heat[_qp] = computeHeatSourceAtTime(x, y, z, _t);
      break;
    case HeatSourceType::LINE:
      _volumetric_heat[_qp] = computeAveragedHeatSource(x, y, z, _t - _dt, _t);
      break;
    case HeatSourceType::MIXED:
      _volumetric_heat[_qp] = computeMixedHeatSource(x, y, z, _t - _dt, _t);
      break;
  }
}

Real
FunctionPathEllipsoidHeatSource::computeHeatSourceAtTime(const Real x,
                                                         const Real y,
                                                         const Real z,
                                                         const Real time)
{
  // The functions that define the path is only time dependent
  const static Point dummy;

  // center of the heat source
  Real x_t = _function_x.value(time, dummy);
  Real y_t = _function_y.value(time, dummy);
  Real z_t = _function_z.value(time, dummy);

  Real dist = std::pow(x - x_t, 2.0) + std::pow(y - y_t, 2.0) + std::pow(z - z_t, 2.0);

  // Gaussian point heat source
  return 2.0 * _P * _eta * _f / libMesh::pi / _r / _r / _r * std::exp(-2.0 * dist / _r / _r);
}

Real
FunctionPathEllipsoidHeatSource::computeAveragedHeatSource(
    const Real x, const Real y, const Real z, const Real time_begin, const Real time_end)
{
  mooseAssert(time_end > time_begin, "Begin time should be smaller than end time.");

  unsigned int num_pts = 5;
  Real Q_integral = 0, Q_integral_old = 0;
  do
  {
    Q_integral_old = Q_integral;
    Real delta_t = (time_end - time_begin) / num_pts;
    Real t0 = time_begin;
    Real Q_begin = computeHeatSourceAtTime(x, y, z, time_begin);
    Real Q_end;
    Q_integral = 0;
    for (unsigned int i = 0; i < num_pts; ++i)
    {
      t0 += delta_t;
      Q_end = computeHeatSourceAtTime(x, y, z, t0);
      // compute integral of Q between t0 and t0 + delta_t
      Q_integral += (Q_begin + Q_end) * delta_t / 2.0;
      // update Q_begin
      Q_begin = Q_end;
    }
    num_pts *= 2;
    // limit to _number_time_integration pts to accelerate the simulation
    if (num_pts > _number_time_integration)
      break;
  } while (!MooseUtils::absoluteFuzzyEqual(Q_integral_old, Q_integral, 1e-4));

  return Q_integral / (time_end - time_begin);
}

Real
FunctionPathEllipsoidHeatSource::computeMixedHeatSource(
    const Real x, const Real y, const Real z, const Real time_begin, const Real time_end)
{
  mooseAssert(time_end > time_begin, "Begin time should be smaller than end time.");

  const static Point dummy;
  // position at time_begin
  Real x_t0 = _function_x.value(time_begin, dummy);
  Real y_t0 = _function_y.value(time_begin, dummy);
  Real z_t0 = _function_z.value(time_begin, dummy);
  Point P_start = Point(x_t0, y_t0, z_t0);
  // position at time_end
  Real x_t = _function_x.value(time_end, dummy);
  Real y_t = _function_y.value(time_end, dummy);
  Real z_t = _function_z.value(time_end, dummy);
  Point P_end = Point(x_t, y_t, z_t);

  // // current position
  // Point P_current = Point(x, y, z);

  // // distances
  // Real SE, SP, PE; // start - current_point - end
  // SE = (P_end - P_start).norm();
  // SP = (P_current - P_start).norm();
  // PE = (P_end - P_current).norm();
  //
  // if (SE < _threshold_length)
  //   return computeHeatSourceAtTime(x, y, z, time_end);
  // else
  // {
  //   if (SP <= _threshold_length || PE <= _threshold_length)
  //     return computeHeatSourceAtTime(x, y, z, time_end);
  //   else
  //     return computeAveragedHeatSource(x, y, z, time_begin, time_end);
  // }

  Real SE = (P_end - P_start).norm();
  if (SE < _threshold_length)
    return computeHeatSourceAtTime(x, y, z, time_end);
  else
    return computeAveragedHeatSource(x, y, z, time_begin, time_end);
}
