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
  params.addClassDescription("Uses Neural Network to control the parameter.");
  params.addRequiredParam<PostprocessorName>(
      "input_pp", "The postprocessor to watch for controlling the specified parameter.");
  params.addRequiredParam<FunctionName>("target",
                                        "The target value 1D time function for the postprocessor");
  params.addParam<std::string>(
      "controled_parameter",
      "The input parameter(s) to control. Specify a single parameter name and all "
      "parameters in all objects matching the name will be updated");
  params.addParam<std::string>("output_pp",
                               "The postprocessor to control. Should be accessed by reference by "
                               "the objects depending on its value.");
  params.addParam<Real>(
      "start_time", -std::numeric_limits<Real>::max(), "The time to start the PID controller at");
  params.addParam<Real>(
      "stop_time", std::numeric_limits<Real>::max(), "The time to stop the PID controller at");

  params.addParam<std::string>(
      "filename", "net.pt", "Filename used to output the neural net parameters.");

  return params;
}

BasicNNControl::BasicNNControl(const InputParameters & parameters)
  : Control(parameters),
    _current(getPostprocessorValueByName(getParam<PostprocessorName>("input_pp"))),
    _parameter_name(getParam<std::string>("controled_parameter")),
    _target(getFunction("target")),
    _start_time(getParam<Real>("start_time")),
    _stop_time(getParam<Real>("stop_time"))
{
  if (isParamValid("controled_parameter") && isParamValid("output_pp"))
    paramError("output_pp",
               "Either a controllable parameter or a postprocessor to control should be specified, "
               "not both.");
  if (!isParamValid("controled_parameter") && !isParamValid("output_pp"))
    mooseError("A controled_parameter or a postprocessor to control should be specified.");
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
    // Compute the delta between the current value of the postprocessor and the desired value
    Real delta = _current - _target.value(_t, dummy);

    // use NN to compute the value

#ifdef TORCH_ENABLED
    auto options = torch::TensorOptions().dtype(at::kDouble);
    std::vector<Real> flattened_data{delta};
    torch::Tensor data_tensor =
        torch::from_blob(flattened_data.data(), {1, 1}, options).to(at::kDouble);
    torch::Tensor prediction = _nn_model->forward(data_tensor);
    auto pred = prediction[0].item<Real>();
    std::cout << "Prediction = " << pred << std::endl;
#endif

    std::cout << "parameter name = " << _parameter_name << std::endl;

    // Set the new value of the controlled parameter
    if (isParamValid("controled_parameter"))
      setControllableValueByName<Real>(_parameter_name, pred);
    else
      _fe_problem.setPostprocessorValueByName(getParam<std::string>("output_pp"), pred);
  }
}
