#include "IncrementGeometryCoupledVarThresholdElementSubdomainModifier.h"

registerMooseObject("SMARTAMApp", IncrementGeometryCoupledVarThresholdElementSubdomainModifier);

InputParameters
IncrementGeometryCoupledVarThresholdElementSubdomainModifier::validParams()
{
  InputParameters params = CoupledVarThresholdElementSubdomainModifier::validParams();

  params.addRequiredParam<std::vector<PostprocessorName>>(
      "postprocessors",
      "The postprocessors which stores x, y, z dimension and the volume of the incremental "
      "geometry size.");

  return params;
}

IncrementGeometryCoupledVarThresholdElementSubdomainModifier::
    IncrementGeometryCoupledVarThresholdElementSubdomainModifier(const InputParameters & parameters)
  : CoupledVarThresholdElementSubdomainModifier(parameters),
    _dim(_fe_problem.mesh().dimension()),
    _postprocessor_names(getParam<std::vector<PostprocessorName>>("postprocessors"))
{
  if (_postprocessor_names.size() != _dim + 1)
    paramError("postprocessors",
               "Please provide the postprocessor names that saves the the x y (z) "
               "dimension and the volume of the "
               "incremental geometry.");
}

void
IncrementGeometryCoupledVarThresholdElementSubdomainModifier::finalize()
{
  CoupledVarThresholdElementSubdomainModifier::finalize();

  // reset the values at every timestep
  _volume = 0.0;
  _max_val.assign(_dim, -DBL_MAX);
  _min_val.assign(_dim, DBL_MAX);

  // obtain the moved elements
  auto elem_range = movedElemsRange();
  for (auto elem : elem_range)
  {
    _volume += elem->volume();
    for (auto node : elem->node_ref_range())
      for (unsigned int i = 0; i < _dim; ++i)
      {
        _min_val[i] = std::min(_min_val[i], node(i));
        _max_val[i] = std::max(_max_val[i], node(i));
      }
  }

  // Get the max/min value across processors
  gatherSum(_volume);
  for (unsigned int i = 0; i < _dim; ++i)
  {
    gatherMin(_min_val[i]);
    gatherMax(_max_val[i]);
  }

  // update postprocessor values
  for (unsigned int i = 0; i < _dim; ++i)
    // if the values are not updated, means no element is moved
    if (MooseUtils::absoluteFuzzyEqual(_max_val[i], DBL_MAX) ||
        MooseUtils::absoluteFuzzyEqual(_min_val[i], DBL_MAX))
      _fe_problem.setPostprocessorValueByName(_postprocessor_names[i], 0.0);
    else
      _fe_problem.setPostprocessorValueByName(_postprocessor_names[i], _max_val[i] - _min_val[i]);

  _fe_problem.setPostprocessorValueByName(_postprocessor_names.back(), _volume);
}
