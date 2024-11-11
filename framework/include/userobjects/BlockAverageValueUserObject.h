//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "ElementUserObject.h"

/**
 * This postprocessor computes a volume integral of the specified
 * variable.
 *
 * Note that specializations of this integral are possible by deriving
 * from this class and overriding computeQpIntegral().
 */
class BlockAverageValueUserObject : public ElementUserObject
{
public:
  static InputParameters validParams();

  BlockAverageValueUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  /// Returns the integral value
  Real averageValue(SubdomainID block) const;

protected:
  virtual Real computeQpIntegral() = 0;
  virtual Real computeIntegral();

  unsigned int _qp;

  // This map will hold the partial sums for each block
  std::map<SubdomainID, Real> _integral_values;

  // This map will hold the partial volume sums for each block
  std::map<SubdomainID, Real> _volume_values;

  // This map will hold our averages for each block
  std::map<SubdomainID, Real> _average_values;
};
