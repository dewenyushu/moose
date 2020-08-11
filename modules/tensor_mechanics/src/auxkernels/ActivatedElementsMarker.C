//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ActivatedElementsMarker.h"

registerMooseObject("TensorMechanicsApp", ActivatedElementsMarker);

template <>
InputParameters
validParams<ActivatedElementsMarker>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Determine activated elements.");
  params.addRequiredCoupledVar("temp_aux",
                               "Temperature aux variable used to determine activated elements.");
  params.addRequiredParam<Real>("melt_temperature", "Melt temperature.");
  params.addParam<bool>("use_avg_temperature",
                        true,
                        "Whether to use the average temperature of an element to decide if an element is activated. The opposite is to activate an element if temperature at any quadrature point is higher than the melt temperature.");

  // params.addParam<UserObjectName>("marker_uo", "Marker UserObject");
  return params;
}

ActivatedElementsMarker::ActivatedElementsMarker(const InputParameters & parameters)
  : AuxKernel(parameters),
    _temp_aux(coupledValue("temp_aux")),
    _melt_temperature(getParam<Real>("melt_temperature")),
    _use_avg_temperature(getParam<bool>("use_avg_temperature"))
// _marker_uo(isParamValid("marker_uo") ? &getUserObjectByName<ActivatedElementsMarkerUO>(
//                                            getParam<UserObjectName>("marker_uo"))
//                                      : nullptr)
{
  // if (_marker_uo)
  //   _marker_map = &(_marker_uo->getActivatedElementsMap());
  // else
  //   _marker_map = nullptr;
}

Real
ActivatedElementsMarker::computeValue()
{
  if (isNodal())
    mooseError("must run on an element variable");

  Real marker_old = _u_old[_qp];
  Real avg_temp = 0;
  bool qp_temp_greater_than_melt = false;
  if (_use_avg_temperature)
  {
    for (unsigned int i = 0; i < _q_point.size(); ++i)
      avg_temp += _temp_aux[i];

    avg_temp /= _q_point.size();
  }
  else
  {
    for (unsigned int i = 0; i < _q_point.size(); ++i)
      if (_temp_aux[i]>_melt_temperature)
      {
        qp_temp_greater_than_melt = true; break;
      }
  }

  if (avg_temp > _melt_temperature || marker_old >= 1 || _current_elem->subdomain_id() != 1 || qp_temp_greater_than_melt)
    return 1;
  else
    return 0;

  // dof_id_type elem_id = _current_elem->id();
  // Real activate_elem = _marker_map->find(elem_id)->second;
  // return activate_elem;
}
