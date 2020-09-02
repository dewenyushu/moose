//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ActivateElementTemp.h"
#include "DisplacedProblem.h"

#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"
#include "libmesh/point.h"

#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/mesh_communication.h"

registerMooseObject("MooseApp", ActivateElementTemp);

template <>
InputParameters
validParams<ActivateElementTemp>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Determine activated elements.");
  params.addRequiredParam<int>("active_subdomain_id", "The active subdomain ID.");
  params.addRequiredParam<std::vector<BoundaryName>>("expand_boundary_name", "The expanded boundary name.");
  params.addParam<FunctionName>("function_x", "The x component heating spot travel path");
  params.addParam<FunctionName>("function_y", "The y component heating spot travel path");
  params.addParam<FunctionName>("function_z", "The z component heating spot travel path");
  params.addParam<bool>(
    "variable_activation",
    false, "Whether to use variable value for element activation. If false, use path activation");
  params.addParam<Real>(
    "activate_value",
    0.0, "The value above which to activate the element");
  params.addCoupledVar(
      "coupled_var",
      "The variable value will be used to decide wether an element whould be activated.");

  return params;
}

ActivateElementTemp::ActivateElementTemp(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _active_subdomain_id(getParam<int>("active_subdomain_id")),
    _expand_boundary_name(getParam<std::vector<BoundaryName>>("expand_boundary_name")),
    _function_x(getFunction("function_x")),
    _function_y(getFunction("function_y")),
    _function_z(getFunction("function_z")),
    _variable_activation(getParam<bool>("variable_activation")),
    _coupled_var(isParamValid("coupled_var") ? & coupledValue("coupled_var"): nullptr),
    _activate_value(getParam<Real>("activate_value"))
{
  if(_variable_activation && _coupled_var == nullptr)
    mooseError("Need a 'coupled_var' defined if variable_activation = true");

  if((!_variable_activation) && _coupled_var!= nullptr)
    mooseWarning("Not using variable activation, so the 'coupled_var' is ignored");

  // add the new boundary and get its boundary id
  _boundary_ids = _mesh.getBoundaryIDs(_expand_boundary_name, true);
}

void
ActivateElementTemp::execute()
{
  /*
    Check if current element is activated
  */

  // activate center (assume position of the activate center is only time dependent)
  Real x_t = _function_x.value(_t, _q_point[0]);
  Real y_t = _function_y.value(_t, _q_point[0]);
  Real z_t = _function_z.value(_t, _q_point[0]);

  Real avg_val = 0.0;
  if(_variable_activation)
  {
      for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
        avg_val +=(* _coupled_var)[qp];
      avg_val /=  _qrule->n_points();
  }

  bool activate_by_var = _variable_activation && avg_val > _activate_value;
  bool activate_by_path = (!_variable_activation) && _current_elem->contains_point( Point (x_t, y_t, z_t));

  if((activate_by_var || activate_by_path) && _current_elem->subdomain_id()!=_active_subdomain_id)
  {
    /*
      _current_elem subdomain id is not assignable
      create a copy of this element from MooseMesh
    */
    dof_id_type ele_id= _current_elem->id();
    Elem * ele = _mesh.elemPtr(ele_id);

    /*
      Add element to the activate subdomain
    */
    ele->subdomain_id()=_active_subdomain_id;
    /*
      Reassign element in the reference mesh while using a displaced mesh
    */
    auto displaced_problem = _fe_problem.getDisplacedProblem();
    if (displaced_problem)
    {
      Elem * disp_ele = displaced_problem->mesh().elemPtr(ele_id);
      disp_ele->subdomain_id()=_active_subdomain_id;

      // Elem * ref_ele = displaced_problem->refMesh().elemPtr(ele_id);
      // ref_ele->subdomain_id()=_active_subdomain_id;
    }
    /*
      Save the newly activated element id for updating boundary info later
    */
    _newly_activated_elem.insert(ele_id);

    // std::cout<<"====>  Neighbor element info:\n";
    // for (auto s : ele->side_index_range())
    // {
    //   if (ele->neighbor_ptr(s))
    //   {
    //     dof_id_type neighbor_ele_id=ele->neighbor_ptr(s)->id();
    //     Elem * neighbor_ele = _mesh.elemPtr(neighbor_ele_id);
    //     neighbor_ele->print_info();
    //   }
    // }
    // std::cout<<"====>  Current element info:\n";
    // ele->print_info();
  }


}

void
ActivateElementTemp::finalize()
{
  /*
    Update boundary info
  */
  updateBoundaryInfo(_mesh);
  /*
    Synchronize ghost element subdomain ID
  */
  SyncSubdomainIds sync(_mesh.getMesh());
  Parallel::sync_dofobject_data_by_id
      (_mesh.getMesh().comm(), _mesh.getMesh().elements_begin(),  _mesh.getMesh().elements_end(), sync);

  auto displaced_problem = _fe_problem.getDisplacedProblem();
  if (displaced_problem)
  {
    updateBoundaryInfo(displaced_problem->mesh());

    SyncSubdomainIds sync_mesh(displaced_problem->mesh().getMesh());
    Parallel::sync_dofobject_data_by_id
        (displaced_problem->mesh().getMesh().comm(), displaced_problem->mesh().getMesh().elements_begin(),  displaced_problem->mesh().getMesh().elements_end(), sync_mesh);
  }


  /*
    Reinit equation systems
  */
  _fe_problem.meshChanged();

  // // print boundary element info
  // ConstBndElemRange & range = *_mesh.getBoundaryElementRange();
  // for (const auto & belem : range)
  // {
  //   std::cout<<"Element ID: "<<belem->_elem->id()
  //             <<"; Boundary name: "<< _mesh.getMesh().get_boundary_info().get_sideset_name(belem->_bnd_id)
  //             <<"; side: "<< belem->_side
  //             <<std::endl;
  // }
}

void ActivateElementTemp::updateBoundaryInfo(MooseMesh & mesh)
{
  // save the removed ghost sides to sync across processors
  std::unordered_map<processor_id_type, std::vector<std::pair<dof_id_type, unsigned int>>>
      ghost_sides_to_remove;

  for (auto ele_id : _newly_activated_elem)
  {
    Elem * ele = mesh.elemPtr(ele_id);
    for (auto s : ele->side_index_range())
    {
      if (ele->neighbor_ptr(s))
      {
        dof_id_type neighbor_ele_id=ele->neighbor_ptr(s)->id();
        Elem * neighbor_ele = mesh.elemPtr(neighbor_ele_id);
        if (neighbor_ele->subdomain_id()!=_active_subdomain_id)
        {
          // add this side to boundary
          mesh.getMesh().get_boundary_info().add_side( ele,  s, _boundary_ids[0]);
        }
        else
        {
          // remove this side from the boundary
          mesh.getMesh().get_boundary_info().remove_side(ele, s,  _boundary_ids[0]);
          // remove the neighbor side from the boundary
          unsigned int neighbor_s = neighbor_ele->which_neighbor_am_i(ele);
          mesh.getMesh().get_boundary_info().remove_side(neighbor_ele, neighbor_s, _boundary_ids[0]);
          if (neighbor_ele->processor_id()!=this->processor_id())
            ghost_sides_to_remove[neighbor_ele->processor_id()].emplace_back(neighbor_ele->id(), neighbor_s);
        }
      }
    }
  }

  // synchronize boundary information across processors
  // mesh.getMesh().get_boundary_info().sync_push_boundary_side_id(ghost_sides_to_remove);
  push_boundary_side_ids(mesh, ghost_sides_to_remove);
  mesh.getMesh().get_boundary_info().sync_pull_boundary_side_id();
  mesh.buildBndElemList();
}

void ActivateElementTemp::push_boundary_side_ids( MooseMesh & mesh
  std::unordered_map<processor_id_type, std::vector<std::pair<dof_id_type, unsigned int>>>
  & elems_to_push)
{
  auto elem_action_functor =
    [& mesh]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, unsigned int>> & received_elem)
    {
      for (const auto & pr : received_elem)
        mesh.getMesh().get_boundary_info().remove_side(mesh.getMesh().elem_ptr(pr.first), pr.second);
    };

  Parallel::push_parallel_vector_data
    (mesh.getMesh().get_boundary_info().comm(), elems_to_push, elem_action_functor);
}
