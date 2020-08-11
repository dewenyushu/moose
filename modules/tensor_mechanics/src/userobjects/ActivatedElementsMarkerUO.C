//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ActivatedElementsMarkerUO.h"
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"

registerMooseObject("MooseApp", ActivatedElementsMarkerUO);

template <>
InputParameters
validParams<ActivatedElementsMarkerUO>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Determine activated elements.");
  params.addRequiredCoupledVar("temp_aux",
                               "Temperature aux variable used to determine activated elements.");
  params.addRequiredParam<Real>("melt_temperature", "Melt temperature.");
  params.addParam<bool>("use_avg_temperature",
                        true,
                        "Whether to use the average temperature of an element to decide if an element is activated. The opposite is to activate an element if temperature at any quadrature point is higher than the melt temperature.");

  return params;
}

ActivatedElementsMarkerUO::ActivatedElementsMarkerUO(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _temp_aux(coupledValue("temp_aux")),
    _melt_temperature(getParam<Real>("melt_temperature")),
    _use_avg_temperature(getParam<bool>("use_avg_temperature"))
{
}

void
ActivatedElementsMarkerUO::execute()
{
  Real marker_old = 0;
  auto marker = _activated_elem_map.find(_current_elem->id());
  if (marker != _activated_elem_map.end())
    marker_old = marker->second;

  Real avg_temp = 0.0;
  bool qp_temp_greater_than_melt = false;
  if (_use_avg_temperature)
  {
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
      avg_temp += _temp_aux[qp];

    avg_temp /= _qrule->n_points();
  }
  else
  {
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
      if (_temp_aux[qp]>_melt_temperature)
      {
        qp_temp_greater_than_melt = true; break;
      }
  }

  if (avg_temp > _melt_temperature || marker_old >= 1.0 || _current_elem->subdomain_id() != 1 || qp_temp_greater_than_melt)
  {
    _activated_elem_map[_current_elem->id()] = 1.0;
    if (marker_old < 1.0)
      _newly_activated_elem.push_back(_current_elem->id());
  }
  else
    _activated_elem_map[_current_elem->id()] = 0.0;
}

void
ActivatedElementsMarkerUO::finalize()
{
  _communicator.set_union(_activated_elem_map);
}
