//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "UpdatedEulerAngle.h"

registerMooseObject("TensorMechanicsApp", UpdatedEulerAngle);

InputParameters
UpdatedEulerAngle::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<bool>(
      "degree_to_radian", true, "Whether to convert euler angles from degree to radian");
  return params;
}

UpdatedEulerAngle::UpdatedEulerAngle(const InputParameters & parameters)
  : Material(parameters),
    _initial_euler_angles(getMaterialPropertyByName<RealVectorValue>("Euler_angles")),
    _update_rotation(getMaterialProperty<RankTwoTensor>("update_rotation")),
    _updated_euler_angle(declareProperty<RealVectorValue>("updated_Euler_angle"))
{
}

void
UpdatedEulerAngle::initQpStatefulProperties()
{
  _updated_euler_angle[_qp] = _initial_euler_angles[_qp];
}

void
UpdatedEulerAngle::computeQpProperties()
{
  computeEulerAngleFromRotationMatrix(_update_rotation[_qp], _updated_euler_angle[_qp]);
}

void
UpdatedEulerAngle::computeEulerAngleFromRotationMatrix(const RankTwoTensor & rot,
                                                       RealVectorValue & euler_angle)
{
  Real phi1, Phi, phi2;

  Phi = std::acos(rot(2, 2));
  if (MooseUtils::absoluteFuzzyEqual(std::abs(rot(2, 2)) - 1.0, 0.0))
  {
    phi1 = 0.0;
    phi2 = std::atan2(rot(0, 1), rot(0, 0));
  }
  else
  {
    Real sPhi = std::sin(Phi);
    phi1 = std::atan2(rot(0, 2) / sPhi, rot(1, 2) / sPhi);
    phi2 = std::atan2(rot(2, 0) / sPhi, -rot(2, 1) / sPhi);
  }

  euler_angle(0) = phi1;
  euler_angle(1) = Phi;
  euler_angle(2) = phi2;

  euler_angle *= 180.0 / pi;
}
