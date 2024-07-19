//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EulerAngleUpdateFromReporter.h"

#include <fstream>

registerMooseObject("SolidMechanicsApp", EulerAngleUpdateFromReporter);

InputParameters
EulerAngleUpdateFromReporter::validParams()
{
  InputParameters params = EulerAngleFileReader::validParams();
  params.addClassDescription("Update Euler angle from reporter value.");
  params.addRequiredParam<ReporterName>(
      "euler_angle_0_name",
      "reporter name for the first component of the Euler angles.  This "
      "uses the reporter syntax <reporter>/<name>.");
  params.addRequiredParam<ReporterName>(
      "euler_angle_1_name",
      "reporter name for the second component of the Euler angles.  This "
      "uses the reporter syntax <reporter>/<name>.");
  params.addRequiredParam<ReporterName>(
      "euler_angle_2_name",
      "reporter name for the third component of the Euler angles.  This "
      "uses the reporter syntax <reporter>/<name>.");
  params.addRequiredParam<ReporterName>("grain_id_name",
                                        "reporter name for the grain IDs.  This "
                                        "uses the reporter syntax <reporter>/<name>.");
  return params;
}

EulerAngleUpdateFromReporter::EulerAngleUpdateFromReporter(const InputParameters & params)
  : EulerAngleFileReader(params),
    _euler_angle_0(
        getReporterValue<std::vector<Real>>("euler_angle_0_name", REPORTER_MODE_REPLICATED)),
    _euler_angle_1(
        getReporterValue<std::vector<Real>>("euler_angle_1_name", REPORTER_MODE_REPLICATED)),
    _euler_angle_2(
        getReporterValue<std::vector<Real>>("euler_angle_2_name", REPORTER_MODE_REPLICATED)),
    _grain_id(getReporterValue<std::vector<Real>>("grain_id_name", REPORTER_MODE_REPLICATED))

{
}

void
EulerAngleUpdateFromReporter::UpdateEulerAngle()
{
  std::cout << "EulerAngleUpdateFromReporter::UpdateEulerAngle()" << std::endl;

  // check sizes of the containers
  if (_grain_id.size() != _euler_angle_0.size() || _grain_id.size() != _euler_angle_1.size() ||
      _grain_id.size() != _euler_angle_2.size())
    paramError("grain_id_name", "Number of reporters' entries do not match.");

  if (_grain_id.size() < _angles.size())
    paramError("grain_id_name",
               "Number of grains from the reporter (",
               _grain_id.size(),
               ") is smaller than the existing number of grains (",
               _angles.size(),
               "). Some existing grains will miss texture information.");

  // zip the grain id with the euler angles
  std::map<int, RealVectorValue> ea_data;
  for (const auto i : index_range(_grain_id))
  {
    int gid = (int)(_grain_id[i]);
    ea_data[gid] = RealVectorValue(_euler_angle_0[i], _euler_angle_1[i], _euler_angle_2[i]);
  }

  // re-assign euler angles based on the new data
  auto max_grain_id = ea_data.rbegin()->first;
  _angles.resize(max_grain_id);

  for (const auto it : ea_data)
  {
    _angles[it.first] = EulerAngles(it.second);
  }
}
