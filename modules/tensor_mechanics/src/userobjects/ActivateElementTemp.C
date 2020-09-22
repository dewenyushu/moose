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
#include "libmesh/dof_map.h"

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
  params.addParam<Real>("activate_distance", 1e-4, "The maximum distance of the activated element to the point on the path.");
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
    _active_subdomain_id(declareRestartableData<int>("active_subdomain_id", getParam<int>("active_subdomain_id"))),
    _inactive_subdomain_id(declareRestartableData<int>("inactive_subdomain_id", getParam<int>("inactive_subdomain_id"))),
    _expand_boundary_name(getParam<std::vector<BoundaryName>>("expand_boundary_name")),
    _function_x(getFunction("function_x")),
    _function_y(getFunction("function_y")),
    _function_z(getFunction("function_z")),
    _activate_distance(declareRestartableData<Real>("activate_distance", getParam<Real>("activate_distance"))),
    _variable_activation(declareRestartableData<bool>("variable_activation", getParam<bool>("variable_activation"))),
    _coupled_var(isParamValid("coupled_var") ? & coupledValue("coupled_var"): nullptr),
    _activate_value(declareRestartableData<Real>("activate_value", getParam<Real>("activate_value")))
{
  if(_variable_activation && _coupled_var == nullptr)
    mooseError("Need a 'coupled_var' defined if variable_activation = true");

  if((!_variable_activation) && _coupled_var!= nullptr)
    mooseWarning("Not using variable activation, so the 'coupled_var' is ignored");

  std::cout<<"active_subdomain_id = "<<_active_subdomain_id<<std::endl;
  std::cout<<"inactive_subdomain_id = "<<_inactive_subdomain_id<<std::endl;
  std::cout<<"expand_boundary_name = "<<_expand_boundary_name[0]<<std::endl;
  std::cout<<"activate_distance = "<<_activate_distance<<std::endl;
  std::cout<<"variable_activation = "<<_variable_activation<<std::endl;
  std::cout<<"activate_value = "<<_activate_value<<std::endl;

  setNewBoundayName();
}

void
ActivateElementTemp::setNewBoundayName()
{
  // add the new boundary and get its boundary id
  _boundary_ids = _mesh.getBoundaryIDs(_expand_boundary_name, true);
  _mesh.setBoundaryName(_boundary_ids[0], _expand_boundary_name[0]);
  _mesh.getMesh().get_boundary_info().sideset_name(_boundary_ids[0]) = _expand_boundary_name[0];
  _mesh.getMesh().get_boundary_info().nodeset_name(_boundary_ids[0]) = _expand_boundary_name[0];

  auto displaced_problem = _fe_problem.getDisplacedProblem();
  if (displaced_problem)
  {
    _disp_boundary_ids=displaced_problem->mesh().getBoundaryIDs(_expand_boundary_name, true);
    displaced_problem->mesh().setBoundaryName(_disp_boundary_ids[0], _expand_boundary_name[0]);
    displaced_problem->mesh().getMesh().get_boundary_info().sideset_name(_disp_boundary_ids[0]) = _expand_boundary_name[0];
    displaced_problem->mesh().getMesh().get_boundary_info().nodeset_name(_disp_boundary_ids[0]) = _expand_boundary_name[0];
  }
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
  bool activate_by_path = (!_variable_activation) && _current_elem->close_to_point( Point (x_t, y_t, z_t), _activate_distance);

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
    Synchronize ghost element subdomain ID
    Note: this needs to be done before updating boundary info because
    updating boundary requires the updated element subdomain ids
  */
  SyncSubdomainIds sync(_mesh.getMesh());
  Parallel::sync_dofobject_data_by_id
      (_mesh.getMesh().comm(), _mesh.getMesh().elements_begin(),  _mesh.getMesh().elements_end(), sync);
  /*
      Update boundary info
  */
  updateBoundaryInfo(_mesh);
  /*
     Similarly for the displaced mesh
  */
  auto displaced_problem = _fe_problem.getDisplacedProblem();
  if (displaced_problem)
  {
    SyncSubdomainIds sync_mesh(displaced_problem->mesh().getMesh());
    Parallel::sync_dofobject_data_by_id
        (displaced_problem->mesh().getMesh().comm(), displaced_problem->mesh().getMesh().elements_begin(),  displaced_problem->mesh().getMesh().elements_end(), sync_mesh);
    updateBoundaryInfo(displaced_problem->mesh());
  }

  /*
    Reinit equation systems
  */
  _fe_problem.meshChanged();

  /*
    Get storage ranges for the newly activated elements and boundary nodes
  */
  ConstElemRange & elem_range = * this->getNewlyActivatedElementRange();
  ConstBndNodeRange & bnd_node_range = * this->getNewlyActivatedBndNodeRange();
  // ConstNodeRange & node_range = * this->getNewlyActivatedNodeRange();
  /*
    Apply initial condition for the newly activated elements
  */
  initSolutions(elem_range, bnd_node_range);

   /*
     Initialize stateful material properties for the newly activated elements
   */
   _fe_problem.initElementStatefulProps(elem_range);

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
      Elem * neighbor_ele = ele->neighbor_ptr(s);
      if (neighbor_ele == nullptr)
      {
        // add this side to boundary
        mesh.getMesh().get_boundary_info().add_side( ele,  s, _boundary_ids[0]);
      }
      else
      {
        if (neighbor_ele->subdomain_id()!=_active_subdomain_id && neighbor_ele->subdomain_id()!=_inactive_subdomain_id)
        {
          // add this side to boundary
          mesh.getMesh().get_boundary_info().add_side( ele,  s, _boundary_ids[0]);
        }
        else
        {
          // remove this side from the boundary
          mesh.getMesh().get_boundary_info().remove_side(ele, s);

          // remove the neighbor side from the boundary
          unsigned int neighbor_s = neighbor_ele->which_neighbor_am_i(ele);
          mesh.getMesh().get_boundary_info().remove_side(neighbor_ele, neighbor_s);

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
    // [&mesh, this]
    [&mesh]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, unsigned int>> & received_elem)
    {
      for (const auto & pr : received_elem)
      {
        // remove the side
        // mesh.getMesh().get_boundary_info().remove_side(mesh.getMesh().elem_ptr(pr.first), pr.second, this->getExpandedBoundaryID());
        mesh.getMesh().get_boundary_info().remove_side(mesh.elemPtr(pr.first), pr.second);
        // // remove the nodes on this side
        // this->remove_bounday_node(mesh, mesh.getMesh().elem_ptr(pr.first)->side_ptr(pr.second));
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

  // std::cout<<"number of newly activated nodes =  "<<_newly_activated_node.size()<<std::endl;

  // create a vector of the newly activated nodes
  std::vector<const BndNode *> nodes;
  std::set<const BndNode *> set_nodes;
  ConstBndNodeRange & bnd_nodes = *_mesh.getBoundaryNodeRange();
  for (auto & bnode : bnd_nodes)
  {
    dof_id_type bnode_id = bnode->_node->id();
    auto it = _newly_activated_node.find(bnode_id);
    if (it != _newly_activated_node.end())
    {
      set_nodes.insert(bnode);
      // std::cout<<"BoundaryNode Id()="<<bnode_id<<std::endl;
    }

  }

  nodes.assign(set_nodes.begin(), set_nodes.end());

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

ConstNodeRange * ActivateElementTemp::getNewlyActivatedNodeRange()
{
  // deletes the object first
  _activated_node_range.reset();

  // create a vector of the newly activated nodes
  std::vector<const Node *> nodes;
  for (auto  elem_id : _newly_activated_elem)
  {
    const Node * const * elem_nodes=_mesh.elemPtr(elem_id)->get_nodes();
    unsigned int n_nodes = _mesh.elemPtr(elem_id)->n_nodes();
    for (unsigned int n =0; n<n_nodes; ++n)
    {
      // check if all the elements connected to this node are newly activated
      const Node * nd = elem_nodes[n];
      if (isNewlyActivated(nd))
        nodes.push_back(nd);
    }
  }

  // Make some fake node iterators defining this vector of
  // nodes
  Node * const * nodepp = const_cast<Node * const * >(nodes.data());
  Node * const * nodeend = nodepp+nodes.size();

  const MeshBase::const_node_iterator nodes_begin =
    MeshBase::const_node_iterator(nodepp,
                                  nodeend,
                                  Predicates::NotNull<Node * const *>());

  const MeshBase::const_node_iterator nodes_end =
    MeshBase::const_node_iterator(nodeend,
                                  nodeend,
                                  Predicates::NotNull<Node * const *>());

  if (!_activated_node_range)
    _activated_node_range = libmesh_make_unique<ConstNodeRange>(
        nodes_begin, nodes_end);

  return _activated_node_range.get();
}

bool ActivateElementTemp::isNewlyActivated(const Node * nd)
{
  const auto & node_to_elem_map = _mesh.nodeToElemMap();
  auto node_to_elem_pair = node_to_elem_map.find(nd->id());
  if (node_to_elem_pair != node_to_elem_map.end())
  {
    const std::vector<dof_id_type> & connected_ele_ids = node_to_elem_pair->second;
    for (auto connected_ele_id : connected_ele_ids)
    {
      // check the connected elements
      if (_mesh.elemPtr(connected_ele_id)->subdomain_id()==_inactive_subdomain_id)
          return false;
      if(_mesh.elemPtr(connected_ele_id)->subdomain_id()==_active_subdomain_id &&
          std::find(_newly_activated_elem.begin(), _newly_activated_elem.end(), connected_ele_id) == _newly_activated_elem.end())
          return false;
    }
  }
  return true;
}

void ActivateElementTemp::remove_bounday_node(MooseMesh & mesh, std::unique_ptr<const Elem> side)
{
  unsigned int n_nodes=side->n_nodes();
  const Node * const * side_nodes = side->get_nodes();
  for (unsigned int n=0; n<n_nodes; n++)
  {
    // only remove side node that does not belong to the newly activated nodes
    if (_newly_activated_node.find (side_nodes[n]->id()) == _newly_activated_node.end())
    {
      // std::cout<<"Remove boundary node id = "<<side_nodes[n]->id()<<std::endl;
      mesh.getMesh().get_boundary_info().remove_node(side_nodes[n], _boundary_ids[0]);
    }
  }
}

void ActivateElementTemp::initSolutions(ConstElemRange & elem_range,
                                        ConstBndNodeRange & bnd_node_range)
{
  // project initial condition to the current solution
  _fe_problem.projectInitialConditionOnCustomRange(elem_range, bnd_node_range);

  NumericVector<Number> & current_solution = * _fe_problem.getNonlinearSystemBase().system().current_local_solution;
  NumericVector<Number> & old_solution = _fe_problem.getNonlinearSystemBase().solutionOld();
  NumericVector<Number> & older_solution = _fe_problem.getNonlinearSystemBase().solutionOlder();

  NumericVector<Number> & current_aux_solution = * _fe_problem.getAuxiliarySystem().system().current_local_solution;
  NumericVector<Number> & old_aux_solution = _fe_problem.getAuxiliarySystem().solutionOld();
  NumericVector<Number> & older_aux_solution = _fe_problem.getAuxiliarySystem().solutionOlder();

  DofMap & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();
  DofMap & dof_map_aux = _fe_problem.getAuxiliarySystem().dofMap();

  // std::cout<<"current sol size = "<<current_solution.size()
  // <<"; old solution size = "<<old_solution.size()
  // <<"; older solution size = "<<older_solution.size()<<std::endl;

  // std::cout<<"current aux sol size = "<<current_aux_solution.size()
  // <<"; old aux solution size = "<<old_aux_solution.size()
  // <<"; older aux solution size = "<<older_aux_solution.size()<<std::endl;

  std::set<dof_id_type> dofs, dofs_aux;
  // get dofs for the newly added elements
  for (auto & elem : elem_range)
  {
    std::vector<dof_id_type> di, di_aux;
    dof_map.dof_indices(elem, di);
    dof_map_aux.dof_indices(elem, di_aux);
    for(unsigned int i =0; i<di.size(); ++i)
      dofs.insert(di[i]);
    for(unsigned int i =0; i<di_aux.size(); ++i)
      dofs_aux.insert(di_aux[i]);

    di.clear();
    di_aux.clear();
  }

  // update solutions
  for (auto dof: dofs)
  {
    old_solution.set(dof, current_solution(dof));
    older_solution.set(dof, current_solution(dof));
  }
  // update aux solutions
  for (auto dof_aux: dofs_aux)
  {
    old_aux_solution.set(dof_aux, current_aux_solution(dof_aux));
    older_aux_solution.set(dof_aux, current_aux_solution(dof_aux));
  }

  dofs.clear();
  dofs_aux.clear();

  current_solution.close();
  old_solution.close();
  older_solution.close();

  current_aux_solution.close();
  old_aux_solution.close();
  older_aux_solution.close();

  _fe_problem.restoreSolutions();
}
