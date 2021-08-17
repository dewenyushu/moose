//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GmshGeneratorBase.h"

class GmshTensionSampleGenerator : public GmshGeneratorBase
{
public:
  static InputParameters validParams();

  GmshTensionSampleGenerator(const InputParameters & parameters);

protected:
  void generateGeometry() override;

  int addCurveCenter(const std::vector<Real> & p1_coord,
                     const std::vector<Real> & p2_coord,
                     const std::vector<Real> & c_coord);

  Real _scale;

  Real _half_bridge_height;
  Real _half_bridge_length;
  Real _gripping_height;
  Real _fillet_radius;
  Real _connection_radius;
  Real _half_sample_length;

  Real _right_circle_x;
  Real _right_circle_y;
  Real _right_circle_r;

  Real _left_circle_x;
  Real _left_circle_y;
  Real _left_circle_r;
};
