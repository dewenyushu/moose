//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"
#include "Function.h"

class ActivateElementTemp;

template <>
InputParameters validParams<ActivateElementTemp>();

class ActivateElementTemp : public ElementUserObject
{
public:
  ActivateElementTemp(const InputParameters & parameters);

  const std::set<dof_id_type> & getNewlyActivatedElements() const
  {
    return _newly_activated_elem;
  };

  BoundaryID getExpandedBoundaryID()
  {
    return _boundary_ids[0];
  }

  void initialize() override{};
  void execute() override;
  void threadJoin(const UserObject & /*uo*/) override{};
  void finalize() override;

protected:
  void updateBoundaryInfo(MooseMesh & mesh);

  void push_boundary_info( MooseMesh & mesh,
    std::unordered_map<processor_id_type, std::vector<std::pair<dof_id_type, unsigned int>>>
    & elems_to_push);

  void remove_bounday_node(MooseMesh & mesh, std::unique_ptr<Elem> side);

  ConstElemRange * getNewlyActivatedElementRange();
  ConstBndNodeRange * getNewlyActivatedBndNodeRange();

  std::set<dof_id_type> _newly_activated_elem;

  /**
   * Ranges for use with threading.
   */
  std::unique_ptr<ConstElemRange> _activated_elem_range;
  std::unique_ptr<ConstBndNodeRange> _activated_bnd_node_range;

  /// activate subdomain ID
  const subdomain_id_type _active_subdomain_id;
  /// expanded boundary name
  const std::vector<BoundaryName> _expand_boundary_name;
  /// expanded boundary IDs
  std::vector<BoundaryID> _boundary_ids;
  /// path of the heat source, x, y, z components
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;
  /// whether to use variable value for element activation
  /// if false, use path activation
  const bool _variable_activation;
  /// varaible to decide wether an element whould be activated
  const VariableValue * _coupled_var;
  /// varaible to decide wether an element whould be activated
  const Real _activate_value;
};
