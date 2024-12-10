//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlockQuaternionMode.h"
#include "MooseMesh.h"

#include "libmesh/mesh_tools.h"

registerMooseObject("SolidMechanicsApp", BlockQuaternionMode);

InputParameters
BlockQuaternionMode::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addParam<unsigned int>("bins", 20, "Number of bins to segregate quaternions");

  return params;
}

BlockQuaternionMode::BlockQuaternionMode(const InputParameters & parameters)
  : ElementUserObject(parameters),
_quaternion_prop(getMaterialProperty<std::vector<Real>>("updated_quaternion")),
_bins(getParam<unsigned int>("bins"))
{
}

Real
BlockQuaternionMode::getBlockValue(SubdomainID block) const
{
  // // Note that we can't use operator[] for a std::map in a const function!
  // if (_average_values.find(block) != _average_values.end())
  //   return _average_values.find(block)->second;

  // mooseError("Unknown block requested for average value!");

  return 0; // To satisfy compilers
}

void
BlockQuaternionMode::initialize()
{
  // // // Explicitly call the initialization routines for our base class
  // // ElementUserObject::initialize();

  // // Set averages to 0 for each block
  // const std::set<SubdomainID> & blocks = _subproblem.mesh().meshSubdomains();

  // for (std::set<SubdomainID>::const_iterator it = blocks.begin(); it != blocks.end(); ++it)
  // {
  //   _integral_values[*it] = 0;
  //   _volume_values[*it] = 0;
  //   _average_values[*it] = 0;
  // }
}

void
BlockQuaternionMode::execute()
{
  // Store quaternion values for the current element
  _quat[_current_elem->subdomain_id()].push_back(_quaternion_prop[0]);
  std::cout<<_current_elem->id()<<"\t"<<_quaternion_prop[0][0]<<"\t"<<_quaternion_prop[0][1]<<"\t"<<_quaternion_prop[0][2]<<"\t"<<_quaternion_prop[0][3]<<std::endl;

  // // // Compute the integral on this element
  // // Real integral_value = computeIntegral();

  // // Add that value to the others we've computed on this subdomain
  // _integral_values[_current_elem->subdomain_id()] += 0.0;

  // // Keep track of the volume of this block
  // _volume_values[_current_elem->subdomain_id()] += _current_elem_volume;
}

void
BlockQuaternionMode::threadJoin(const UserObject & y)
{
  // // ElementUserObject::threadJoin(y);

  // // We are joining with another class like this one so do a cast so we can get to it's data
  // const BlockQuaternionMode & bav = dynamic_cast<const BlockQuaternionMode &>(y);

  // for (std::map<SubdomainID, Real>::const_iterator it = bav._integral_values.begin();
  //      it != bav._integral_values.end();
  //      ++it)
  //   _integral_values[it->first] += it->second;

  // for (std::map<SubdomainID, Real>::const_iterator it = bav._volume_values.begin();
  //      it != bav._volume_values.end();
  //      ++it)
  //   _volume_values[it->first] += it->second;

  // for (std::map<SubdomainID, Real>::const_iterator it = bav._average_values.begin();
  //      it != bav._average_values.end();
  //      ++it)
  //   _average_values[it->first] += it->second;
}

void
BlockQuaternionMode::finalize()
{
  // // Loop over the integral values and sum them up over the processors
  // for (std::map<SubdomainID, Real>::iterator it = _integral_values.begin();
  //      it != _integral_values.end();
  //      ++it)
  //   gatherSum(it->second);

  // // Loop over the volumes and sum them up over the processors
  // for (std::map<SubdomainID, Real>::iterator it = _volume_values.begin();
  //      it != _volume_values.end();
  //      ++it)
  //   gatherSum(it->second);

  // // Now everyone has the correct data so everyone can compute the averages properly:
  // for (std::map<SubdomainID, Real>::iterator it = _average_values.begin();
  //      it != _average_values.end();
  //      ++it)
  // {
  //   SubdomainID id = it->first;
  //   if (MooseUtils::absoluteFuzzyEqual(_volume_values[id], 0.0))
  //     _average_values[id] = 0.0;
  //   else
  //     _average_values[id] = _integral_values[id] / _volume_values[id];
  // }

  std::cout<<_quat.size()<<std::endl;
}
