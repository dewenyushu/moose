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
    _volumetric_heat(declareADProperty<Real>("volumetric_heat"))
{
}

void
FunctionPathEllipsoidHeatSource::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  // The functions that define the path is only time dependent
  const static Point dummy;

  // center of the heat source
  Real x_t = _function_x.value(_t, dummy);
  Real y_t = _function_y.value(_t, dummy);
  Real z_t = _function_z.value(_t, dummy);

  Real dist = std::pow(x - x_t, 2.0) + std::pow(y - y_t, 2.0) + std::pow(z - z_t, 2.0);

  _volumetric_heat[_qp] =
      2.0 * _P * _eta * _f / libMesh::pi / _r / _r * std::exp(-2.0 * dist / _r / _r);
}
