//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GaussianHeatSourceBase.h"


/**
 * Double ellipsoid heat source distribution.
 */
class VelocityGaussianHeatSource : public GaussianHeatSourceBase
{
public:
  static InputParameters validParams();

  VelocityGaussianHeatSource(const InputParameters & parameters);

protected:
  virtual void
  computeHeatSourceCenterAtTime(Real & x, Real & y, Real & z, const Real & time) override;

  virtual void computeHeatSourceMovingSpeedAtTime(const Real & time) override;

  /// previous time
  Real _prev_time;

  /// position at previous time
  Real _x_prev, _y_prev, _z_prev;

  /// function of scanning speed along three directions
  const Function & _function_vx;
  const Function & _function_vy;
  const Function & _function_vz;

  /// stores the velocity at the current step
  Real _vx, _vy, _vz;
};
