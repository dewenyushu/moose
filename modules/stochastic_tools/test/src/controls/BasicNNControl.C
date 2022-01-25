//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BasicNNControl.h"
#include "Function.h"
#include "Transient.h"

registerMooseObject("MooseApp", BasicNNControl);

InputParameters
BasicNNControl::validParams()
{
  InputParameters params = Control::validParams();
  params.addClassDescription(
      "Sets the value of a 'Real' input parameter (or postprocessor) based on a Proportional "
      "Integral Derivative control of a postprocessor to match a target a target value.");
  params.addRequiredParam<PostprocessorName>(
      "postprocessor", "The postprocessor to watch for controlling the specified parameter.");
  params.addRequiredParam<FunctionName>("target",
                                        "The target value 1D time function for the postprocessor");
  params.addRequiredParam<Real>("K_integral", "The coefficient multiplying the integral term");
  params.addRequiredParam<Real>("K_proportional",
                                "The coefficient multiplying the difference term");
  params.addRequiredParam<Real>("K_derivative", "The coefficient multiplying the derivative term");
  params.addParam<std::string>(
      "parameter",
      "The input parameter(s) to control. Specify a single parameter name and all "
      "parameters in all objects matching the name will be updated");
  params.addParam<std::string>("parameter_pp",
                               "The postprocessor to control. Should be accessed by reference by "
                               "the objects depending on its value.");
  params.addParam<Real>(
      "start_time", -std::numeric_limits<Real>::max(), "The time to start the PID controller at");
  params.addParam<Real>(
      "stop_time", std::numeric_limits<Real>::max(), "The time to stop the PID controller at");
  params.addParam<bool>(
      "reset_every_timestep",
      false,
      "Reset the PID integral when changing timestep, for coupling iterations within a timestep");
  params.addParam<bool>("reset_integral_windup",
                        true,
                        "Reset the PID integral when the error crosses zero and the integral is "
                        "larger than the error.");

  params.addParam<std::string>(
      "filename", "net.pt", "Filename used to output the neural net parameters.");

  return params;
}

BasicNNControl::BasicNNControl(const InputParameters & parameters)
  : Control(parameters),
    _current(getPostprocessorValueByName(getParam<PostprocessorName>("postprocessor"))),
    _target(getFunction("target")),
    _Kint(getParam<Real>("K_integral")),
    _Kpro(getParam<Real>("K_proportional")),
    _Kder(getParam<Real>("K_derivative")),
    _start_time(getParam<Real>("start_time")),
    _stop_time(getParam<Real>("stop_time")),
    _reset_every_timestep(getParam<bool>("reset_every_timestep")),
    _reset_integral_windup(getParam<bool>("reset_integral_windup")),
    _integral(0),
    _integral_old(0),
    _value_old(0),
    _t_step_old(-1),
    _old_delta(0)
// #ifdef TORCH_ENABLED
//     ,
//     _nn_model(std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet>)
// #endif

{
  if (isParamValid("parameter") && isParamValid("parameter_pp"))
    paramError("parameter_pp",
               "Either a controllable parameter or a postprocessor to control should be specified, "
               "not both.");
  if (!isParamValid("parameter") && !isParamValid("parameter_pp"))
    mooseError("A parameter or a postprocessor to control should be specified.");
  if (isParamValid("parameter") && _reset_every_timestep)
    paramError(
        "reset_every_timestep",
        "Resetting the PID every time step is only supported using controlled postprocessors");
  if (!dynamic_cast<Transient *>(_app.getExecutioner()))
    mooseError("BasicNNControl is only supported by a Transient Executioner. If using a Steady "
               "Executioner, add a minimum number of Picard iterations and comment this error.");
#ifdef TORCH_ENABLED
  std::vector<unsigned int> nurons_per_layer{5, 5};
  // Initialize neural network model
  _nn_model = std::make_shared<StochasticTools::LibtorchSimpleNeuralNet>(
      "nn_model.pt", 1, 2, nurons_per_layer, 1);

#endif
}

void
BasicNNControl::execute()
{
  Point dummy;

  if (_t >= _start_time && _t < _stop_time)
  {
    // Get the current value of the controllable parameter
    Real value;
    if (isParamValid("parameter"))
      value = getControllableValue<Real>("parameter");
    else
      value = getPostprocessorValueByName(getParam<std::string>("parameter_pp"));

    // Save integral and controlled value at each time step
    // if the solver fails, a smaller time step will be used but _t_step is unchanged
    if (_t_step != _t_step_old)
    {
      // Reset the error integral if PID is only used within each timestep
      if (_reset_every_timestep)
        _integral = 0;

      _integral_old = _integral;
      _value_old = value;
      _t_step_old = _t_step;
      _old_delta = 0;
    }

    // If there were coupling/Picard iterations during the transient and they failed,
    // we need to reset the controlled value and the error integral to their initial value at the
    // beginning of the coupling process
    if (_app.getExecutioner()->picardSolve().numPicardIts() == 1)
    {
      _integral = _integral_old;
      value = _value_old;
    }

    // Compute the delta between the current value of the postprocessor and the desired value
    Real delta = _current - _target.value(_t, dummy);

    ////-- use NN to compute the value, START

    // If desired, reset integral of the error if the error crosses zero
    if (_reset_integral_windup && delta * _old_delta < 0)
      _integral = 0;

    // Compute the three error terms and add them to the controlled value
    _integral += delta * _dt;
    value += _Kint * _integral + _Kpro * delta;
    if (_dt > 0)
      value += _Kder * delta / _dt;

    // Set the new value of the postprocessor
    if (isParamValid("parameter"))
      setControllableValue<Real>("parameter", value);
    else
      _fe_problem.setPostprocessorValueByName(getParam<std::string>("parameter_pp"), value);

#ifdef TORCH_ENABLED
    auto options = torch::TensorOptions().dtype(at::kDouble);
    std::vector<Real> flattened_data{delta};
    torch::Tensor data_tensor =
        torch::from_blob(flattened_data.data(), {1, 1}, options).to(at::kDouble);
    torch::Tensor prediction = _nn_model->forward(data_tensor);
    auto pred = prediction[0].item<Real>();
    std::cout << "Prediction = " << pred << std::endl;
#endif

    /// -- use NN to compute the value, END

    // Keep track of the previous delta for integral windup control
    _old_delta = delta;
  }
}
