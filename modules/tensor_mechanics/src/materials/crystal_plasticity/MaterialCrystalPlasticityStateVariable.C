//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MaterialCrystalPlasticityStateVariable.h"

#include <fstream>

registerMooseObject("TensorMechanicsApp", MaterialCrystalPlasticityStateVariable);

InputParameters
MaterialCrystalPlasticityStateVariable::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<FileName>(
      "state_variable_file_name",
      "",
      "Name of the file containing the initial values of slip system resistances");
  MooseEnum intvar_read_options("file_input inline_input user_input", "inline_input");
  params.addParam<MooseEnum>(
      "intvar_read_type",
      intvar_read_options,
      "Read from options for initial value of internal variables: Default from .i file");
  params.addParam<std::vector<unsigned int>>("groups",
                                             "To group the initial values on different "
                                             "slip systems 'format: [start end)', i.e.'0 "
                                             "4 8 11' groups 0-3, 4-7 and 8-11 ");
  params.addParam<std::vector<Real>>("group_values",
                                     "The initial values corresponding to each "
                                     "group, i.e. '0.0 1.0 2.0' means 0-2 = 0.0, "
                                     "4-7 = 1.0 and 8-11 = 2.0 ");
  params.addParam<std::vector<std::string>>("uo_state_var_evol_rate_comp_name",
                                            "Name of state variable evolution rate component "
                                            "property: Same as state variable evolution rate "
                                            "component user object specified in input file.");
  params.addParam<Real>("zero", 0.0, "Numerical zero for interval variable");
  params.addParam<std::vector<Real>>("scale_factor", "Scale factor of individual component.");

  params.addRequiredParam<unsigned int>("variable_size", "The variable's size.");
  params.addClassDescription(
      "Crystal plasticity state variable class.  Override the virtual functions in your class");
  return params;
}

MaterialCrystalPlasticityStateVariable::MaterialCrystalPlasticityStateVariable(
    const InputParameters & parameters)
  : Material(parameters),
    _num_mat_state_var_evol_rate_comps(
        parameters.get<std::vector<std::string>>("uo_state_var_evol_rate_comp_name").size()),
    _mat_prop_state_var(getMaterialProperty<std::vector<Real>>(_name)),
    _state_variable_file_name(getParam<FileName>("state_variable_file_name")),
    _intvar_read_type(getParam<MooseEnum>("intvar_read_type")),
    _groups(getParam<std::vector<unsigned int>>("groups")),
    _group_values(getParam<std::vector<Real>>("group_values")),
    _zero(getParam<Real>("zero")),
    _scale_factor(getParam<std::vector<Real>>("scale_factor")),
    // Stateful material
    _variable_size(getParam<unsigned int>("variable_size")),
    _values(declareProperty<std::vector<Real>>("property_name")),
    _values_old(getMaterialPropertyOld<std::vector<Real>>("property_name"))
{
  if (_scale_factor.size() != _num_mat_state_var_evol_rate_comps)
    mooseError(
        "MaterialCrystalPlasticityStateVariable: Scale factor should be have the same size of "
        "evolution rate components.");

  _mat_prop_state_var_evol_rate_comps.resize(_num_mat_state_var_evol_rate_comps);

  for (unsigned int i = 0; i < _num_mat_state_var_evol_rate_comps; ++i)
    _mat_prop_state_var_evol_rate_comps[i] = &getMaterialProperty<std::vector<Real>>(
        parameters.get<std::vector<std::string>>("uo_state_var_evol_rate_comp_name")[i]);
}

void
MaterialCrystalPlasticityStateVariable::initQpStatefulProperties()
{
  if (_groups.size() <= 0)
    mooseError(
        "MaterialCrystalPlasticityStateVariable: Error in reading initial state variable values: "
        "Specify input in .i file or in state_variable file");
  else if (_groups.size() != (_group_values.size() + 1))
    mooseError("MaterialCrystalPlasticityStateVariable: The size of the groups and group_values "
               "does not match.");

  for (unsigned int i = 0; i < _groups.size() - 1; ++i)
  {
    unsigned int is, ie;

    is = _groups[i];
    ie = _groups[i + 1] - 1;

    if (is > ie)
      mooseError("MaterialCrystalPlasticityStateVariable: Start index is = ",
                 is,
                 " should be greater than end index ie = ",
                 ie,
                 " in state variable read");

    for (unsigned int j = is; j <= ie; ++j)
      _values[_qp][j] = _group_values[i];
  }
}

void
MaterialCrystalPlasticityStateVariable::readInitialValueFromFile() const
{
  MooseUtils::checkFileReadable(_state_variable_file_name);

  std::ifstream file;
  file.open(_state_variable_file_name.c_str());

  for (unsigned int i = 0; i < _variable_size; ++i)
    if (!(file >> _values[_qp][i]))
      mooseError(
          "Error MaterialCrystalPlasticityStateVariable: Premature end of state_variable file");

  file.close();
}

void
MaterialCrystalPlasticityStateVariable::readInitialValueFromInline() const
{
}

void
MaterialCrystalPlasticityStateVariable::provideInitialValueByUser() const
{
  mooseError("Error MaterialCrystalPlasticityStateVariable: User has to overwrite "
             "'provideInitialValueByUser' function"
             "in order to provide specific initial values based on quadrature point location.");
}

void
MaterialCrystalPlasticityStateVariable::computeQpProperties()
{
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    _values[_qp][i] = 0.0;
    for (unsigned int j = 0; j < _num_mat_state_var_evol_rate_comps; j++)
      _values[_qp][i] += (*_mat_prop_state_var_evol_rate_comps[j])[_qp][i] * _dt * _scale_factor[j];
  }

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    if (_values_old[_qp][i] < _zero && _values[_qp][i] < 0.0)
      _values[_qp][i] = _values_old[_qp][i];
    else
      _values[_qp][i] = _values_old[_qp][i] + _values[_qp][i];

    // if (_values[i] < 0.0)
    //   return false;
  }
  // return true;
}
