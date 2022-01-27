//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#ifdef TORCH_ENABLED
#include <torch/torch.h>
#include "LibtorchSimpleNeuralNet.h"
#endif

#include "Control.h"

/**
 * A time-dependent control of an input parameter or a postprocessor, which aims at
 * making a postprocessor match a desired value.
 */
class BasicNNControl : public Control, public torch::nn::Module
{
public:
  static InputParameters validParams();

  /**
   * Class constructor
   * @param parameters Input parameters for this Control object
   */
  BasicNNControl(const InputParameters & parameters);

  virtual void execute() override;

private:
  /// The current value of the target postprocessor
  const PostprocessorValue & _current;
  /// Name of the controlled parameter
  const std::string & _parameter_name;
  /// The target 1D time-dependent function for the postprocessor
  const Function & _target;
  /// The time to start the controller on
  const Real _start_time;
  /// The time to stop using the controller on
  const Real _stop_time;

#ifdef TORCH_ENABLED
  std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet> _nn_model;
#endif
};
