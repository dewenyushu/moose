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
#include "Function.h"
#include "MooseEnum.h"

/**
 * Double ellipsoid heat source distribution.
 */
class FunctionPathEllipsoidHeatSource : public Material
{
public:
  static InputParameters validParams();

  FunctionPathEllipsoidHeatSource(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  Real computeHeatSourceAtTime(const Real x, const Real y, const Real z, const Real time);

  Real computeAveragedHeatSource(
      const Real x, const Real y, const Real z, const Real time_begin, const Real time_end);

  Real computeMixedHeatSource(
      const Real x, const Real y, const Real z, const Real time_begin, const Real time_end);

  /// power
  const Real _P;
  /// process efficienty
  const Real _eta;
  /// effective radius
  const Real _r;
  /// scaling factor
  const Real _f;
  /// path of the heat source, x, y, z components
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;

  /// type of heat source
  const enum class HeatSourceType { POINT, LINE, MIXED } _heat_source_type;

  const Real _threshold_length;

  ADMaterialProperty<Real> & _volumetric_heat;
};
