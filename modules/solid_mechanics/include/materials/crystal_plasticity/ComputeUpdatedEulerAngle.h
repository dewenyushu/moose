//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "EulerAngles.h"

class ComputeUpdatedEulerAngle : public Material
{
public:
  static InputParameters validParams();

  ComputeUpdatedEulerAngle(const InputParameters & parameters);

  void initQpStatefulProperties() override;
  void computeQpProperties() override;

private:
  void computePropertiesFromRotationMatrix(const RankTwoTensor & rot,
                                           RealVectorValue & euler_angle,
                                           std::vector<Real> & quaternion);

  // updated rotation tensor
  const MaterialProperty<RankTwoTensor> & _updated_rotation;
  // updated euler angle that is computed in thie class
  MaterialProperty<RealVectorValue> & _updated_euler_angle;
  // updated quaternion
  MaterialProperty<std::vector<Real>> & _updated_quaternion;
};
