//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EulerAngleMagnitude.h"

registerMooseObject("TensorMechanicsApp", EulerAngleMagnitude);

InputParameters
EulerAngleMagnitude::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<bool>(
      "degree_to_radian", true, "Whether to convert euler angles from degree to radian");
  return params;
}

EulerAngleMagnitude::EulerAngleMagnitude(const InputParameters & parameters)
  : Material(parameters),
    _degree_to_radian(getParam<bool>("degree_to_radian")),
    // _initial_euler_angles(getMaterialPropertyByName<RealVectorValue>("Euler_angles")),
    _update_rotation(getMaterialProperty<RankTwoTensor>("update_rot")),
    _updated_euler_angle(declareProperty<RealVectorValue>("updated_Euler_angle"))
{
}

void
EulerAngleMagnitude::computeQpProperties()
{
  computeEulerAngleFromRotationMatrix(_update_rotation[_qp], _updated_euler_angle[_qp]);

  if (!_degree_to_radian)
    _updated_euler_angle[_qp] *= 180.0;
}

void
EulerAngleMagnitude::computeEulerAngleFromRotationMatrix(const RankTwoTensor & rot,
                                                         RealVectorValue & euler_angle)
{
  euler_angle.zero();
  if (MooseUtils::absoluteFuzzyEqual(rot(2, 0) - 1.0, 0.0) ||
      MooseUtils::absoluteFuzzyEqual(rot(2, 0) + 1.0, 0.0))
  {
    euler_angle(2) = 0;
    if (MooseUtils::absoluteFuzzyEqual(rot(2, 0) + 1.0, 0.0))
    {
      euler_angle(1) = pi / 2.0;
      euler_angle(0) = euler_angle(2) + std::atan2(rot(0, 1), rot(0, 2));
    }
    else
    {
      euler_angle(1) = -pi / 2.0;
      euler_angle(0) = -euler_angle(2) + std::atan2(-rot(0, 1), -rot(0, 2));
    }
  }
  else
  {
    // we have two sets of possible euler angles in this case, use only one of them
    euler_angle(1) = -std::asin(rot(2, 0));
    euler_angle(0) =
        std::atan2(rot(2, 1) / std::cos(euler_angle(1)), rot(2, 2) / std::cos(euler_angle(1)));
    euler_angle(2) =
        std::atan2(rot(1, 0) / std::cos(euler_angle(1)), rot(0, 0) / std::cos(euler_angle(1)));
  }
}
