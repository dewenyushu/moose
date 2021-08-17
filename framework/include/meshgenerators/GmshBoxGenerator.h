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

class GmshBoxGenerator : public GmshGeneratorBase
{
public:
  static InputParameters validParams();

  GmshBoxGenerator(const InputParameters & parameters);

protected:
  void generateGeometry() override;

  Real _x, _y, _z;
};
