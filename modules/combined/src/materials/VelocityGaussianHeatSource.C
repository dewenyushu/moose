//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VelocityGaussianHeatSource.h"

registerMooseObject("CombinedApp", VelocityGaussianHeatSource);

InputParameters
VelocityGaussianHeatSource::validParams()
{
  InputParameters params = GaussianHeatSourceBase::validParams();
  params.addParam<Real>("x0", 0, "The x component of the initial center of the heating spot");
  params.addParam<Real>("y0", 0, "The y component of the initial center of the heating spot");
  params.addParam<Real>("z0", 0, "The z component of the initial center of the heating spot");

  params.addParam<FunctionName>("function_vx", 0, "The function of x component of the center of the heating spot moving speed");
  params.addParam<FunctionName>("function_vy", 0, "The function of y component of the center of the heating spot moving speed");
  params.addParam<FunctionName>("function_vz", 0, "The function of z component of the center of the heating spot moving speed");

  params.addClassDescription("Double ellipsoid volumetric source heat with moving velocity.");

  return params;
}

VelocityGaussianHeatSource::VelocityGaussianHeatSource(const InputParameters & parameters)
  : GaussianHeatSourceBase(parameters),
    _prev_time(0),
    _x_prev(getParam<Real>("x0")),
    _y_prev(getParam<Real>("y0")),
    _z_prev(getParam<Real>("z0")),
    _function_vx(getFunction("function_vx")),
    _function_vy(getFunction("function_vy")),
    _function_vz(getFunction("function_vz"))
{
}

void
VelocityGaussianHeatSource::computeHeatSourceCenterAtTime(Real & x,
                                                          Real & y,
                                                          Real & z,
                                                          const Real & time)
{
  Real delta_t;
  if (time > _prev_time)
  {
    const static Point dummy;
    _vx = _function_vx.value(time, dummy);
    _vy = _function_vy.value(time, dummy);
    _vz = _function_vz.value(time, dummy);

    delta_t = time - _prev_time;

    x = _x_prev + _vx * delta_t;
    y = _y_prev + _vy * delta_t;
    z = _z_prev + _vz * delta_t;

    // march forward time
    _prev_time = time;
    // update position
    _x_prev += _vx * delta_t;
    _y_prev += _vy * delta_t;
    _z_prev += _vz * delta_t;
  }
  else
  {
    x = _x_prev;
    y = _y_prev;
    z = _z_prev;
  }
}

void
VelocityGaussianHeatSource::computeHeatSourceMovingSpeedAtTime(const Real & /*time*/)
{
  _scan_speed = std::sqrt(_vx * _vx + _vy * _vy + _vz * _vz);
}
