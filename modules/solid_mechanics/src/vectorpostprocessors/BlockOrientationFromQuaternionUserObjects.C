//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlockOrientationFromQuaternionUserObjects.h"
#include "MooseMesh.h"
#include "Assembly.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "libmesh/quadrature.h"
#include "EulerAngles.h"

registerMooseObject("SolidMechanicsApp", BlockOrientationFromQuaternionUserObjects);

InputParameters
BlockOrientationFromQuaternionUserObjects::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();

  params.addRequiredParam<std::vector<UserObjectName>>("quaternion_average_uos", "List of BlockAverage user objects for quaternion components.");

  params.addParam<bool>(
      "degree_to_radian", false, "Whether to convert euler angles from degree to radian.");

  params.addClassDescription("Output the Euler angle for each block computed from average of quaternions.");
  return params;
}

BlockOrientationFromQuaternionUserObjects::BlockOrientationFromQuaternionUserObjects(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _mesh(_subproblem.mesh()),
    _uo_names(getParam<std::vector<UserObjectName>>("quaternion_average_uos")),
    _num_cols(4), // add one colum for the subdomain ID
    _num_rows(_mesh.meshSubdomains().size())
{
  _output_vector.resize(_num_cols);
  _uos.resize(4); // 4 quaternion components
  for (const auto j : make_range(_num_cols))
  {
    if (j == 0)
      _output_vector[j] = &declareVector("subdomain_id");
    else
    {
      _output_vector[j] = &declareVector("avg_ea"+std::to_string(j)); // can change
    }
  }

  for (int j=0; j<4; ++j)
        _uos[j] = &getUserObjectByName<BlockAverage>(_uo_names[j]);
}

void
BlockOrientationFromQuaternionUserObjects::initialize()
{
  for (const auto j : make_range(_num_cols))
  {
    _output_vector[j]->clear();
    _output_vector[j]->resize(_num_rows, 0.0);
  }
}

void
BlockOrientationFromQuaternionUserObjects::finalize()
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
BlockOrientationFromQuaternionUserObjects::execute()
{
  int row = 0;
  for (const auto sid : _mesh.meshSubdomains())
  {
    // get Euler angle for each subdomain
    Eigen::Quaternion<Real> q(_uos[0]->averageValue(sid),_uos[1]->averageValue(sid),_uos[2]->averageValue(sid),_uos[3]->averageValue(sid));
    // construct Euler angle from Quaternion
    EulerAngles ea(q);
    // convert EulerAngles to RealVectorValue
    RealVectorValue euler_angle = (RealVectorValue)ea;

    if (getParam<bool>("degree_to_radian"))
      euler_angle *= pi / 180.0;

    for (const auto col : make_range(_num_cols))
    {
      if (col == 0)
        (*_output_vector[col])[row] = sid;
      else
        (*_output_vector[col])[row] = euler_angle(col-1);
    }
    row++;
  }
}
