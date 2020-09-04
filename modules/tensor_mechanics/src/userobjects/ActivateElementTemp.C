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
  params.addParam<int>("inactive_subdomain_id", -1 , "The inactivate subdomain ID, i.e., the subdomain that you want to keep the same.");
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
    _inactive_subdomain_id(getParam<int>("inactive_subdomain_id")),
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


  // initializeBoundary(_mesh);
  // /*
  //   Synchronize ghost element subdomain ID
  // */
  // SyncSubdomainIds sync(_mesh.getMesh());
  // Parallel::sync_dofobject_data_by_id
  //     (_mesh.getMesh().comm(), _mesh.getMesh().elements_begin(),  _mesh.getMesh().elements_end(), sync);
  //
  // auto displaced_problem = _fe_problem.getDisplacedProblem();
  // if (displaced_problem)
  // {
  //   initializeBoundary(displaced_problem->mesh());
  //
  //   SyncSubdomainIds sync_mesh(displaced_problem->mesh().getMesh());
  //   Parallel::sync_dofobject_data_by_id
  //       (displaced_problem->mesh().getMesh().comm(), displaced_problem->mesh().getMesh().elements_begin(),  displaced_problem->mesh().getMesh().elements_end(), sync_mesh);
  // }
  //
  // /*
  //   Notify the mesh about the change
  // */
  // _mesh.meshChanged();
  // displaced_problem->mesh().meshChanged();
}

void
ActivateElementTemp::initializeBoundary(MooseMesh & mesh)
{
  // We want the active subdomain to be empty, but it is often not the case.
  // So we want to check the initial boundary info and update the moving boundary

  // save the removed ghost sides and associated nodes to sync across processors
  std::unordered_map<processor_id_type, std::vector<std::pair<dof_id_type, unsigned int>>>
      ghost_sides_to_remove;

  for (auto & ele: *mesh.getActiveLocalElementRange())
  {
    if (ele->subdomain_id() == _active_subdomain_id)
    {
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
            // add the nodes on this side to the _newly_activated_node
            unsigned int n_nodes=ele->side_ptr(s)->n_nodes();
            for (unsigned int n =0; n<n_nodes; ++n)
              _newly_activated_node.insert(ele->side_ptr(s)->node_ptr(n)->id());
          }
          else
          {
            // remove this side from the boundary
            mesh.getMesh().get_boundary_info().remove_side(ele, s,  _boundary_ids[0]);
            remove_bounday_node(mesh, ele->side_ptr(s));

            // remove the neighbor side from the boundary
            unsigned int neighbor_s = neighbor_ele->which_neighbor_am_i(ele);
            mesh.getMesh().get_boundary_info().remove_side(neighbor_ele, neighbor_s, _boundary_ids[0]);
            remove_bounday_node(mesh, neighbor_ele->side_ptr(neighbor_s));

            if (neighbor_ele->processor_id()!=this->processor_id())
              ghost_sides_to_remove[neighbor_ele->processor_id()].emplace_back(neighbor_ele->id(), neighbor_s);
          }
        }
      }
    }
  }

  // synchronize boundary information across processors
  push_boundary_info(mesh, ghost_sides_to_remove);
  mesh.getMesh().get_boundary_info().parallel_sync_side_ids();
  mesh.getMesh().get_boundary_info().parallel_sync_node_ids();
  // mesh.buildBndElemList();
  // mesh.buildNodeList();
  mesh.update();

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

  if((activate_by_var || activate_by_path) && _current_elem->subdomain_id()!=_active_subdomain_id
      &&  _current_elem->subdomain_id()!=_inactive_subdomain_id )
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

    }
    /*
      Save the newly activated element id for updating boundary info later
    */
    _newly_activated_elem.insert(ele_id);

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

  /*
    Initialize stateful material properties for the newly activated elements
  */
  ConstElemRange & elem_range = * this->getNewlyActivatedElementRange();
  _fe_problem.initElementStatefulProps(elem_range);
  /*
    Apply initial condition for the newly activated elements and nodes
  */
  _fe_problem.projectInitialConditionOnElemRange(elem_range);

  ConstBndNodeRange & bnd_node_range = * this->getNewlyActivatedBndNodeRange();
  _fe_problem.projectInitialConditionOnNodeRange(bnd_node_range);

  /*
    Clear the list
  */
  _newly_activated_elem.clear();
  _newly_activated_node.clear();
}

void ActivateElementTemp::updateBoundaryInfo(MooseMesh & mesh)
{
  // save the removed ghost sides and associated nodes to sync across processors
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
          // add the nodes on this side to the _newly_activated_node
          unsigned int n_nodes=ele->side_ptr(s)->n_nodes();
          for (unsigned int n =0; n<n_nodes; ++n)
            _newly_activated_node.insert(ele->side_ptr(s)->node_ptr(n)->id());
        }
        else
        {
          // remove this side from the boundary
          mesh.getMesh().get_boundary_info().remove_side(ele, s,  _boundary_ids[0]);
          remove_bounday_node(mesh, ele->side_ptr(s));

          // remove the neighbor side from the boundary
          unsigned int neighbor_s = neighbor_ele->which_neighbor_am_i(ele);
          mesh.getMesh().get_boundary_info().remove_side(neighbor_ele, neighbor_s, _boundary_ids[0]);
          remove_bounday_node(mesh, neighbor_ele->side_ptr(neighbor_s));

          if (neighbor_ele->processor_id()!=this->processor_id())
            ghost_sides_to_remove[neighbor_ele->processor_id()].emplace_back(neighbor_ele->id(), neighbor_s);
        }
      }
    }
  }

  // synchronize boundary information across processors
  push_boundary_info(mesh, ghost_sides_to_remove);
  mesh.getMesh().get_boundary_info().parallel_sync_side_ids();
  mesh.getMesh().get_boundary_info().parallel_sync_node_ids();
  // mesh.buildBndElemList();
  // mesh.buildNodeList();
  mesh.update();
}

void ActivateElementTemp::push_boundary_info( MooseMesh & mesh,
  std::unordered_map<processor_id_type, std::vector<std::pair<dof_id_type, unsigned int>>>
  & elems_to_push)
{
  auto elem_action_functor =
    [&mesh, this]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, unsigned int>> & received_elem)
    {
      for (const auto & pr : received_elem)
      {
        // remove the side
        mesh.getMesh().get_boundary_info().remove_side(mesh.getMesh().elem_ptr(pr.first), pr.second, this->getExpandedBoundaryID());
        // remove the nodes on this side
        this->remove_bounday_node(mesh, mesh.getMesh().elem_ptr(pr.first)->side_ptr(pr.second));
      }
    };

  Parallel::push_parallel_vector_data
    (mesh.getMesh().get_boundary_info().comm(), elems_to_push, elem_action_functor);
}


ConstElemRange * ActivateElementTemp::getNewlyActivatedElementRange()
{
  // deletes the object first
  _activated_elem_range.reset();

  // create a vector of the newly activated elements
  std::vector<Elem*> elems;
  for (auto elem_id : _newly_activated_elem)
    elems.push_back(_mesh.elemPtr(elem_id));

  // Make some fake element iterators defining this vector of
  // elements
  Elem * const * elempp = const_cast<Elem * const * >(elems.data());
  Elem * const * elemend = elempp+elems.size();

  const MeshBase::const_element_iterator elems_begin =
    MeshBase::const_element_iterator(elempp,
                                     elemend,
                                     Predicates::NotNull<Elem * const *>());

  const MeshBase::const_element_iterator elems_end =
    MeshBase::const_element_iterator(elemend,
                                     elemend,
                                     Predicates::NotNull<Elem * const *>());
  if (!_activated_elem_range)
    _activated_elem_range = libmesh_make_unique<ConstElemRange>(
        elems_begin, elems_end);

  return _activated_elem_range.get();
}

ConstBndNodeRange * ActivateElementTemp::getNewlyActivatedBndNodeRange()
{
  // deletes the object first
  _activated_bnd_node_range.reset();

  // create a vector of the newly activated nodes
  std::vector<const BndNode *> nodes;
  ConstBndNodeRange & bnd_nodes = *_mesh.getBoundaryNodeRange();
  for (auto & bnode : bnd_nodes)
  {
    dof_id_type bnode_id = bnode->_node->id();
    auto it = _newly_activated_node.find(bnode_id);
    if (it != _newly_activated_node.end())
    {
      nodes.push_back(bnode);
      // std::cout<<"BoundaryNode Id()="<<bnode_id<<std::endl;
    }

  }

  // Make some fake element iterators defining this vector of
  // nodes
  BndNode * const * nodepp = const_cast<BndNode * const * >(nodes.data());
  BndNode * const * nodeend = nodepp+nodes.size();

  const MooseMesh::const_bnd_node_iterator nodes_begin =
    MooseMesh::const_bnd_node_iterator(nodepp,
                                  nodeend,
                                  Predicates::NotNull<BndNode * const *>());

  const MooseMesh::const_bnd_node_iterator nodes_end =
    MooseMesh::const_bnd_node_iterator(nodeend,
                                  nodeend,
                                  Predicates::NotNull<BndNode * const *>());

  if (!_activated_bnd_node_range)
    _activated_bnd_node_range = libmesh_make_unique<ConstBndNodeRange>(
        nodes_begin, nodes_end);

  return _activated_bnd_node_range.get();
}


void ActivateElementTemp::remove_bounday_node(MooseMesh & mesh, std::unique_ptr<const Elem> side)
{
  unsigned int n_nodes=side->n_nodes();
  const Node * const * side_nodes = side->get_nodes();
  for (unsigned int n=0; n<n_nodes; n++)
    mesh.getMesh().get_boundary_info().remove_node(side_nodes[n], _boundary_ids[0]);
}
