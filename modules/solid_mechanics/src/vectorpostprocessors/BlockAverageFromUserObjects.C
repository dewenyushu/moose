//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlockAverageFromUserObjects.h"
#include "MooseMesh.h"
#include "Assembly.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "libmesh/quadrature.h"

registerMooseObject("SolidMechanicsApp", BlockAverageFromUserObjects);

InputParameters
BlockAverageFromUserObjects::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();

  params.addRequiredParam<std::vector<UserObjectName>>("block_average_uos", "List of BlockAverage user objects.");

  params.addClassDescription("Output the block average of variables provided by BlockAverage user objects.");
  return params;
}

BlockAverageFromUserObjects::BlockAverageFromUserObjects(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _mesh(_subproblem.mesh()),
    _uo_names(getParam<std::vector<UserObjectName>>("block_average_uos")),
    _num_cols(_uo_names.size() + 1), // add one colum for the subdomain ID
    _num_rows(_mesh.meshSubdomains().size())
{
  _output_vector.resize(_num_cols);
  _uos.resize(_num_cols-1);
  for (const auto j : make_range(_num_cols))
  {
    if (j == 0)
      _output_vector[j] = &declareVector("subdomain_id");
    else
    {
      _output_vector[j] = &declareVector(_uo_names[j - 1]);
      _uos[j-1] = &getUserObjectByName<BlockAverage>(_uo_names[j - 1]);
    }

  }
}

void
BlockAverageFromUserObjects::initialize()
{
  for (const auto j : make_range(_num_cols))
  {
    _output_vector[j]->clear();
    _output_vector[j]->resize(_num_rows, 0.0);
  }
}

void
BlockAverageFromUserObjects::finalize()
{
  // parallel communication
  for (const auto row : make_range(_num_rows))
  {
    for (const auto col : make_range(_num_cols))
    {
        _communicator.max((*_output_vector[col])[row]);
    }
  }
}

void
BlockAverageFromUserObjects::execute()
{
  int row = 0;
  for (const auto sid : _mesh.meshSubdomains())
  {
    for (const auto col : make_range(_num_cols))
    {
      if (col == 0)
        (*_output_vector[col])[row] = sid;
      else
        (*_output_vector[col])[row] = _uos[col-1]->averageValue(sid);
    }
    row++;
  }
}
