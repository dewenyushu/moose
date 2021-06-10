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

class UpdatedEulerAngle : public Material
{
public:
  static InputParameters validParams();

  UpdatedEulerAngle(const InputParameters & parameters);

  void initQpStatefulProperties() override;
  void computeQpProperties() override;

private:
  void computeEulerAngleFromRotationMatrix(const RankTwoTensor & rot, RealVectorValue & euler_angle);

  const MaterialProperty<RealVectorValue> & _initial_euler_angles;

  // updated rotation tensor
  const MaterialProperty<RankTwoTensor> & _update_rotation;
  MaterialProperty<RealVectorValue> & _updated_euler_angle;
};
