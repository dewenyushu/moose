#pragma once

#include <cfloat>

#include "CoupledVarThresholdElementSubdomainModifier.h"
#include "libmesh/dense_vector.h"

class IncrementGeometryCoupledVarThresholdElementSubdomainModifier
  : public CoupledVarThresholdElementSubdomainModifier
{
public:
  static InputParameters validParams();

  IncrementGeometryCoupledVarThresholdElementSubdomainModifier(const InputParameters & parameters);

protected:
  virtual void finalize() override;

private:
  // Problem dimension
  const unsigned int _dim;

  // Minimum and maximum coordinate values of the incremental geometry. This is updated at every
  // timestep.
  std::vector<Real> _max_val, _min_val;

  // Volume of the changes elements
  Real _volume;

  // Postprocessors that saves and outputs the dimension values calculated from _max_val and
  // _min_val
  std::vector<PostprocessorName> _postprocessor_names;
};
