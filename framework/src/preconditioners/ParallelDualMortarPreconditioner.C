//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ParallelDualMortarPreconditioner.h"

// MOOSE includes
#include "FEProblem.h"
#include "MooseUtils.h"
#include "MooseVariableFE.h"
#include "NonlinearSystem.h"
#include "ComputeJacobianBlocksThread.h"
#include "MooseEnum.h"

#include "libmesh/coupling_matrix.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_base.h"
#include "libmesh/variable.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/parallel_object.h"
#include "libmesh/boundary_info.h"

registerMooseObjectAliased("MooseApp", ParallelDualMortarPreconditioner, "PDMP");

defineLegacyParams(ParallelDualMortarPreconditioner);

InputParameters
ParallelDualMortarPreconditioner::validParams()
{
  InputParameters params = MoosePreconditioner::validParams();

  params.addClassDescription("Dual mortar preconditioner (PDMP) condenses out the Lagrange "
                             "multipliers from the Jacobian matrix "
                             "and recover a system with only the primal unkowns in order to allow "
                             "for a broader range of solvers/preconditioners.");

  params.addParam<std::vector<NonlinearVariableName>>(
      "off_diag_row",
      "The off diagonal row you want to add into the matrix, it will be associated "
      "with an off diagonal column from the same position in off_diag_colum.");
  params.addParam<std::vector<NonlinearVariableName>>(
      "off_diag_column",
      "The off diagonal column you want to add into the matrix, it will be "
      "associated with an off diagonal row from the same position in "
      "off_diag_row.");
  params.addParam<std::vector<NonlinearVariableName>>(
      "coupled_groups",
      "List multiple space separated groups of comma separated variables. "
      "Off-diagonal jacobians will be generated for all pairs within a group.");
  params.addParam<bool>("full",
                        false,
                        "Set to true if you want the full set of couplings.  Simply "
                        "for convenience so you don't have to set every "
                        "off_diag_row and off_diag_column combination.");
  params.addRequiredParam<BoundaryID>("primary_boundary", "Primary side of the contact interface.");
  params.addRequiredParam<BoundaryID>("secondary_boundary",
                                      "Secondary side of the contact interface.");
  params.addRequiredParam<SubdomainID>("primary_subdomain", "Primary subdomain.");
  params.addRequiredParam<SubdomainID>("secondary_subdomain", "Secondary subdomain.");
  params.addRequiredParam<std::vector<std::string>>("preconditioner", "Preconditioner type.");

  return params;
}

ParallelDualMortarPreconditioner::ParallelDualMortarPreconditioner(const InputParameters & params)
  : MoosePreconditioner(params),
    Preconditioner<Number>(MoosePreconditioner::_communicator),
    _nl(_fe_problem.getNonlinearSystemBase()),
    _mesh(&_fe_problem.mesh()),
    _dofmap(&_nl.system().get_dof_map()),
    _n_vars(_nl.nVariables()),
    _K2ci(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _K2cc(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _D(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _M(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _MDinv(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _J_condensed(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _primary_boundary(getParam<BoundaryID>("primary_boundary")),
    _secondary_boundary(getParam<BoundaryID>("secondary_boundary")),
    _primary_subdomain(getParam<SubdomainID>("primary_subdomain")),
    _secondary_subdomain(getParam<SubdomainID>("secondary_subdomain")),
    _save_dofs(false),
    _init_timer(registerTimedSection("init", 2)),
    _apply_timer(registerTimedSection("apply", 1))
{
  // check if SubdomainID & BoundaryID are valid
  if (_mesh->meshSubdomains().find(_secondary_subdomain) == _mesh->meshSubdomains().end())
    mooseError("secondary subdomain ID ", _secondary_subdomain, " does not exist.");
  if (_mesh->meshSubdomains().find(_primary_subdomain) == _mesh->meshSubdomains().end())
    mooseError("primary subdomain ID ", _primary_subdomain, " does not exist.");
  if (_mesh->getBoundaryIDs().find(_secondary_boundary) == _mesh->getBoundaryIDs().end())
    mooseError("Secondary boundary ID ", _secondary_boundary, " does not exist.");
  if (_mesh->getBoundaryIDs().find(_primary_boundary) == _mesh->getBoundaryIDs().end())
    mooseError("Secondary boundary ID ", _primary_boundary, " does not exist.");

  // PC type
  const std::vector<std::string> & pc_type = getParam<std::vector<std::string>>("preconditioner");
  if (pc_type.size() > 1)
    mooseWarning("We only use one preconditioner type in PDMP, the ",
                 pc_type[0],
                 " preconditioner is utilized.");
  _pre_type = Utility::string_to_enum<PreconditionerType>(pc_type[0]);

  std::unique_ptr<CouplingMatrix> cm = libmesh_make_unique<CouplingMatrix>(_n_vars);
  bool full = getParam<bool>("full");

  if (!full)
  {
    // put 1s on diagonal
    for (unsigned int i = 0; i < _n_vars; i++)
      (*cm)(i, i) = 1;

    // off-diagonal entries from the off_diag_row and off_diag_column parameters
    std::vector<std::vector<unsigned int>> off_diag(_n_vars);
    for (unsigned int i = 0;
         i < getParam<std::vector<NonlinearVariableName>>("off_diag_row").size();
         i++)
    {
      unsigned int row =
          _nl.getVariable(0, getParam<std::vector<NonlinearVariableName>>("off_diag_row")[i])
              .number();
      unsigned int column =
          _nl.getVariable(0, getParam<std::vector<NonlinearVariableName>>("off_diag_column")[i])
              .number();
      (*cm)(row, column) = 1;
    }

    // off-diagonal entries from the coupled_groups parameters
    std::vector<NonlinearVariableName> groups =
        getParam<std::vector<NonlinearVariableName>>("coupled_groups");
    for (unsigned int i = 0; i < groups.size(); ++i)
    {
      std::vector<NonlinearVariableName> vars;
      MooseUtils::tokenize<NonlinearVariableName>(groups[i], vars, 1, ",");
      for (unsigned int j = 0; j < vars.size(); ++j)
        for (unsigned int k = j + 1; k < vars.size(); ++k)
        {
          unsigned int row = _nl.getVariable(0, vars[j]).number();
          unsigned int column = _nl.getVariable(0, vars[k]).number();
          (*cm)(row, column) = 1;
          (*cm)(column, row) = 1;
        }
    }
  }
  else
  {
    for (unsigned int i = 0; i < _n_vars; i++)
      for (unsigned int j = 0; j < _n_vars; j++)
        (*cm)(i, j) = 1;
  }

  _fe_problem.setCouplingMatrix(std::move(cm));

  _nl.attachPreconditioner(this);
}

ParallelDualMortarPreconditioner::~ParallelDualMortarPreconditioner() { this->clear(); }

void
ParallelDualMortarPreconditioner::getDofVarSubdomain()
{
  _dof_sets.resize(_n_vars);
  for (unsigned int vn = 0; vn < _n_vars; vn++)
  {
    Variable var = _nl.system().variable(vn);
    // loop through active subdomains of this variable
    for (auto it = var.active_subdomains().begin(); it != var.active_subdomains().end(); ++it)
    {
      SubdomainID sd = (SubdomainID)(*it);
      std::set<dof_id_type> dofs;
      // check dofs of each element
      ConstElemRange * active_elems = _mesh->getActiveElementRange();
      for (const auto & elem : *active_elems)
      {
        if (elem->subdomain_id() == sd)
        {
          std::vector<dof_id_type> di;
          _dofmap->dof_indices(elem, di, vn);
          dofs.insert(di.begin(), di.end());
        }
      }
      _dof_sets[vn].insert(make_pair(sd, dofs));
    }
  }
}

void
ParallelDualMortarPreconditioner::getLocalDofVarSubdomain()
{
  _local_dof_sets.resize(_n_vars);
  for (unsigned int vn = 0; vn < _n_vars; vn++)
  {
    Variable var = _nl.system().variable(vn);
    // loop through active subdomains of this variable
    for (auto it = var.active_subdomains().begin(); it != var.active_subdomains().end(); ++it)
    {
      SubdomainID sd = (SubdomainID)(*it);
      std::set<dof_id_type> dofs;
      // check dofs of each element
      ConstElemRange * active_local_elems = _mesh->getActiveLocalElementRange();
      for (const auto & elem : *active_local_elems)
      {
        if (elem->subdomain_id() == sd)
        {
          std::vector<dof_id_type> di;
          _dofmap->dof_indices(elem, di, vn);
          for (auto index : di)
          {
            if (_dofmap->local_index(index))
              dofs.insert(index);
          }
        }
      }
      _local_dof_sets[vn].insert(make_pair(sd, dofs));
    }
  }
}

void
ParallelDualMortarPreconditioner::getDofVarInterior()
{
  getDofVarSubdomain();
  _dof_sets_interior.resize(_n_vars);
  for (unsigned int vn = 0; vn < _n_vars; vn++)
  {
    Variable var = _nl.system().variable(vn);
    // loop through active subdomains of this variable
    for (auto it = var.active_subdomains().begin(); it != var.active_subdomains().end(); ++it)
    {
      SubdomainID sd = (SubdomainID)(*it);
      std::set<dof_id_type> dofs = _dof_sets[vn][sd];
      // remove dofs on the interface
      for (auto it_primary = _dof_sets_primary[vn].begin();
           it_primary != _dof_sets_primary[vn].end();
           ++it_primary)
        dofs.erase(*it_primary);
      for (auto it_secondary = _dof_sets_secondary[vn].begin();
           it_secondary != _dof_sets_secondary[vn].end();
           ++it_secondary)
        dofs.erase(*it_secondary);

      std::vector<dof_id_type> vec_dofs(dofs.begin(), dofs.end());
      _dof_sets_interior[vn].insert(make_pair(sd, vec_dofs));
    }
  }
}

void
ParallelDualMortarPreconditioner::getLocalDofVarInterior()
{
  getLocalDofVarSubdomain();
  _local_dof_sets_interior.resize(_n_vars);
  for (unsigned int vn = 0; vn < _n_vars; vn++)
  {
    Variable var = _nl.system().variable(vn);
    // loop through active subdomains of this variable
    for (auto it = var.active_subdomains().begin(); it != var.active_subdomains().end(); ++it)
    {
      SubdomainID sd = (SubdomainID)(*it);
      std::set<dof_id_type> dofs = _local_dof_sets[vn][sd];
      // remove dofs on the interface
      for (auto it_primary = _local_dof_sets_primary[vn].begin();
           it_primary != _local_dof_sets_primary[vn].end();
           ++it_primary)
        dofs.erase(*it_primary);
      for (auto it_secondary = _local_dof_sets_secondary[vn].begin();
           it_secondary != _local_dof_sets_secondary[vn].end();
           ++it_secondary)
        dofs.erase(*it_secondary);

      std::vector<dof_id_type> vec_dofs(dofs.begin(), dofs.end());
      _local_dof_sets_interior[vn].insert(make_pair(sd, vec_dofs));
    }
  }
}

void
ParallelDualMortarPreconditioner::getDofVarInterface()
{
  _dof_sets_primary.resize(_n_vars);
  _dof_sets_secondary.resize(_n_vars);
  _dof_sets_secondary_unsorted.resize(_n_vars);
  for (unsigned int vn = 0; vn < _n_vars; vn++)
  {
    // loop over boundary nodes
    ConstBndNodeRange & range = *_mesh->getBoundaryNodeRange();
    std::vector<dof_id_type> di;
    for (const auto & bnode : range)
    {
      const Node * node_bdry = bnode->_node;
      BoundaryID boundary_id = bnode->_bnd_id;

      if (boundary_id == _secondary_boundary)
      {
        _dofmap->dof_indices(node_bdry, di, vn);
        for (auto it = di.begin(); it != di.end(); ++it)
          if (std::find(_dof_sets_secondary[vn].begin(), _dof_sets_secondary[vn].end(), *it) ==
              _dof_sets_secondary[vn].end())
            _dof_sets_secondary[vn].push_back(*it);
      }

      if (boundary_id == _primary_boundary)
      {
        _dofmap->dof_indices(node_bdry, di, vn);
        for (auto it = di.begin(); it != di.end(); ++it)
          if (std::find(_dof_sets_primary[vn].begin(), _dof_sets_primary[vn].end(), *it) ==
              _dof_sets_primary[vn].end())
            _dof_sets_primary[vn].push_back(*it);
      }
    }
    // save an unsorted copy of the indices
    _dof_sets_secondary_unsorted[vn].resize(_dof_sets_secondary[vn].size());
    _dof_sets_secondary_unsorted[vn] = _dof_sets_secondary[vn];
    // sort indices
    std::sort(_dof_sets_primary[vn].begin(), _dof_sets_primary[vn].end());
    // std::sort(_dof_sets_secondary[vn].begin(), _dof_sets_secondary[vn].end());
  }
  // save a map for getting dof of u2c with lm that are associated with the same node
  for (auto i : index_range(_dof_sets_secondary[1]))
    _lm_to_u2c.insert(std::make_pair(_dof_sets_secondary[1][i], _dof_sets_secondary[0][i]));
  // sort the lm indices
  std::sort(_dof_sets_secondary[1].begin(), _dof_sets_secondary[1].end());
  // save corresponding u2c indices
  for (auto i : index_range(_dof_sets_secondary[0]))
    _dof_sets_secondary[0][i] = _lm_to_u2c[_dof_sets_secondary[1][i]];

#ifdef DEBUG
  std::cout << "gu2c = ";
  for (auto i : _dof_sets_secondary[0])
    std::cout << i << ", ";
  std::cout << std::endl;

  std::cout << "glm = ";
  for (auto i : _dof_sets_secondary[1])
    std::cout << i << ", ";
  std::cout << std::endl;
#endif
}

void
ParallelDualMortarPreconditioner::getLocalDofVarInterface()
{
  _local_dof_sets_primary.resize(_n_vars);
  _local_dof_sets_secondary.resize(_n_vars);
  _local_dof_sets_secondary_unsorted.resize(_n_vars);
  for (unsigned int vn = 0; vn < _n_vars; vn++)
  {
    // loop over boundary nodes
    ConstBndNodeRange & range = *_mesh->getBoundaryNodeRange();
    std::vector<dof_id_type> di;
    for (const auto & bnode : range)
    {
      const Node * node_bdry = bnode->_node;
      BoundaryID boundary_id = bnode->_bnd_id;

      if (boundary_id == _secondary_boundary)
      {
        _dofmap->dof_indices(node_bdry, di, vn);
        for (auto it = di.begin(); it != di.end(); ++it)
          if (std::find(_local_dof_sets_secondary[vn].begin(),
                        _local_dof_sets_secondary[vn].end(),
                        *it) == _local_dof_sets_secondary[vn].end() &&
              _dofmap->local_index(*it))
            _local_dof_sets_secondary[vn].push_back(*it);
      }

      if (boundary_id == _primary_boundary)
      {
        _dofmap->dof_indices(node_bdry, di, vn);
        for (auto it = di.begin(); it != di.end(); ++it)
          if (std::find(_local_dof_sets_primary[vn].begin(),
                        _local_dof_sets_primary[vn].end(),
                        *it) == _local_dof_sets_primary[vn].end() &&
              _dofmap->local_index(*it))
            _local_dof_sets_primary[vn].push_back(*it);
      }
    }
    // save a unsorted copy for the lower_D dofs
    _local_dof_sets_secondary_unsorted[vn].resize(_local_dof_sets_secondary[vn].size());
    _local_dof_sets_secondary_unsorted[vn] = _local_dof_sets_secondary[vn];
    // sort indices for primary
    std::sort(_local_dof_sets_primary[vn].begin(), _local_dof_sets_primary[vn].end());
    // secondary dofs for lm needs have a one-one correspondance of u
    // so it needs to be sorted separately
  }
  if (_local_dof_sets_secondary[0].size() != _local_dof_sets_secondary[1].size())
    mooseWarning("lm and u2c do not have same number of local dofs");

  std::sort(_local_dof_sets_secondary[1].begin(), _local_dof_sets_secondary[1].end());
  // save corresponding u2c indices
  for (auto i : index_range(_local_dof_sets_secondary[0]))
    _local_dof_sets_secondary[0][i] = _lm_to_u2c[_local_dof_sets_secondary[1][i]];

#ifdef DEBUG
  std::cout << "u2c = ";
  for (auto i : _local_dof_sets_secondary[0])
    std::cout << i << ", ";
  std::cout << std::endl;

  std::cout << "lm = ";
  for (auto i : _local_dof_sets_secondary[1])
    std::cout << i << ", ";
  std::cout << std::endl;
#endif
}

void
ParallelDualMortarPreconditioner::condenseSystem()
{
  std::vector<dof_id_type> u1c = _local_dof_sets_primary[0];
  std::vector<dof_id_type> u2c = _local_dof_sets_secondary[0];

  std::vector<dof_id_type> lm = _local_dof_sets_secondary[1];

  std::vector<dof_id_type> u1i = _local_dof_sets_interior[0][_primary_subdomain];
  std::vector<dof_id_type> u2i = _local_dof_sets_interior[0][_secondary_subdomain];

  // global indices
  std::vector<dof_id_type> gu1c = _dof_sets_primary[0];
  std::vector<dof_id_type> gu2c = _dof_sets_secondary[0];

  std::vector<dof_id_type> glm = _dof_sets_secondary[1];

  std::vector<dof_id_type> gu1i = _dof_sets_interior[0][_primary_subdomain];
  std::vector<dof_id_type> gu2i = _dof_sets_interior[0][_secondary_subdomain];

#ifdef DEBUG
  std::cout << "u1c = ";
  for (auto i : u1c)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "lm = ";
  for (auto i : lm)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "glm = ";
  for (auto i : glm)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "u2c = ";
  for (auto i : u2c)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "gu2c = ";
  for (auto i : gu2c)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "_rows = ";
  for (auto i : _rows)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "_cols = ";
  for (auto i : _cols)
    std::cout << i << " ";
  std::cout << std::endl;
#endif

  _matrix->create_submatrix(*_M, u1c, lm);

#ifdef DEBUG
  // std::cout << "_M = \n";
  // _M->print_personal();
  std::cout << "Norms of _M_transpose: l1-norm = " << _M->l1_norm()
            << "; infinity-norm = " << _M->linfty_norm() << "\n";
  // << "; frobenius_norm = " << _M->frobenius_norm() << std::endl;
#endif

  _matrix->create_submatrix(*_K2ci, u2c, u2i);
  _matrix->create_submatrix(*_K2cc, u2c, u2c);

#ifdef DEBUG
  std::cout << "Norms of _K2ci: l1-norm = " << _K2ci->l1_norm()
            << "; infinity-norm = " << _K2ci->linfty_norm() << "\n";
  // << "; frobenius_norm = " << _K2ci->frobenius_norm() << std::endl;
  std::cout << "Norms of _K2cc: l1-norm = " << _K2cc->l1_norm()
            << "; infinity-norm = " << _K2cc->linfty_norm() << "\n";
  // << "; frobenius_norm = " << _K2cc->frobenius_norm() << std::endl;
#endif

  // invert _D:
  // _D should be strictly diagonal if dual_mortar approach is utilized
  // so we only need to compute the reciprocal number of the diagonal entries
  // to save memory, no new matrix is created
  _matrix->create_submatrix(*_D, u2c, lm);
#ifdef DEBUG
  std::cout << "Norms of _D: l1-norm = " << _D->l1_norm()
            << "; infinity-norm = " << _D->linfty_norm() << "\n";
  // << "; frobenius_norm = " << _D->frobenius_norm() << std::endl;
#endif
  auto diag_D = NumericVector<Number>::build(MoosePreconditioner::_communicator);
  // Allocate storage
  diag_D->init(_D->n(), _D->local_n(), false, PARALLEL);
  // Fill entries
  for (numeric_index_type i = _D->row_start(); i < _D->row_stop(); ++i)
    diag_D->set(_u2c_to_idx[gu2c[i]], (*_D)(i, _u2c_to_idx[gu2c[i]]));

  _D->zero();

  for (numeric_index_type i = _D->row_start(); i < _D->row_stop(); ++i)
    if (!MooseUtils::absoluteFuzzyEqual((*diag_D)(i), 0.0))
      _D->set(i, _u2c_to_idx[gu2c[i]], 1.0 / (*diag_D)(i));
  _D->close();

#ifdef DEBUG
  // std::cout << "_Dinv = \n";
  // _D->print_personal();
#endif

  // compute MDinv=_M*_D
  _M->matrix_matrix_mult(*_D, *_MDinv); // (should use empty initializer for _MDinv)

#ifdef DEBUG
  // std::cout << "_MDinv =\n";
  // _MDinv->print_personal();
  std::cout << "Norms of _MDinv: l1-norm = " << _MDinv->l1_norm()
            << "; infinity-norm = " << _MDinv->linfty_norm() << "\n";
  // << "; frobenius_norm = " << _MDinv->frobenius_norm() << std::endl;
#endif

  // initialize _J_condensed
  _J_condensed->init(_grows.size(), _gcols.size(), _rows.size(), _cols.size());
  MatSetOption(_J_condensed->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  _matrix->create_submatrix_nosort(*_J_condensed, _grows, _gcols);
  // _matrix->create_submatrix(*_J_condensed, _rows, _cols);

#ifdef DEBUG
  std::cout << "Norms of _J_condensed after cureate_submatrix: l1-norm = "
            << _J_condensed->l1_norm() << "; infinity-norm = " << _J_condensed->linfty_norm()<<"\n";
            // << "; frobenius_norm = " << _J_condensed->frobenius_norm() << std::endl;
  // _J_condensed->print_personal();
#endif

  // compute changed parts: MDinv*K2ci, MDinv*K2cc
  std::unique_ptr<PetscMatrix<Number>> MDinvK2ci(
      libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
      MDinvK2cc(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator));

  // _matrix->create_submatrix(
  //     *MDinvK2ci, u1c, u2i); // get MDinvK2ci initialized (should use empty initializer here)
  //
  // _matrix->create_submatrix(
  //     *MDinvK2cc, u1c, u2c); // get MDinvK2cc initialized (should use empty initializer here)
  _MDinv->matrix_matrix_mult(*_K2ci, *MDinvK2ci);
  _MDinv->matrix_matrix_mult(*_K2cc, *MDinvK2cc);

  MDinvK2ci->close();
  MDinvK2cc->close();

#ifdef DEBUG
  // std::cout << "Submatrix MDinvK2ci =\n";
  // MDinvK2ci->print_personal();
  std::cout << "Norms of MDinvK2ci: l1-norm = " << MDinvK2ci->l1_norm()
            << "; infinity-norm = " << MDinvK2ci->linfty_norm() << "\n";
  // << "; Frobenius-norm = " << MDinvK2ci->frobenius_norm() << std::endl;
  // std::cout << "Submatrix MDinvK2cc =\n";
  // MDinvK2cc->print_personal();
  std::cout << "Norms of MDinvK2cc: l1-norm = " << MDinvK2cc->l1_norm()
            << "; infinity-norm = " << MDinvK2cc->linfty_norm() << "\n";
  // << "; Frobenius-norm = " << MDinvK2cc->frobenius_norm() << std::endl;
#endif

  // add changed parts to _J_condensed
  // original system row_id: u1c
  // original system col_id: u2i, u2c
  std::vector<numeric_index_type> row_id_cond, col_id_cond_u2i, col_id_cond_u2c;
  std::map<numeric_index_type, numeric_index_type> row_id_cond_mp, row_id_cond_u2i_mp,
      row_id_cond_u2c_mp, col_id_cond_u2i_mp, col_id_cond_u2c_mp;

#ifdef DEBUG
  std::cout << "grows = ";
  for (auto i : _grows)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "gcols = ";
  for (auto i : _gcols)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "gu2c = ";
  for (auto i : gu2c)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "gu2i = ";
  for (auto i : gu2i)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "gu1c = ";
  for (auto i : gu1c)
    std::cout << i << " ";
  std::cout << std::endl;
#endif

  for (auto it : index_range(gu1c))
  {
    numeric_index_type lid = static_cast<numeric_index_type>(it);
    auto it_row = find(_grows.begin(), _grows.end(), gu1c[it]);
    if (it_row != _grows.end())
    {
      numeric_index_type gid = std::distance(_grows.begin(), it_row);
      row_id_cond_mp.insert(std::make_pair(lid, gid));
      if (lid >= MDinvK2ci->row_start() && lid < MDinvK2ci->row_stop())
        row_id_cond_u2i_mp.insert(std::make_pair(lid, gid));
      if (lid >= MDinvK2cc->row_start() && lid < MDinvK2cc->row_stop())
        row_id_cond_u2c_mp.insert(std::make_pair(lid, gid));
    }
    else
      mooseError("DOF ", gu1c[it], " does not exist in the rows of the condensed system");
  }

  // global cols
  for (auto it : index_range(gu2i))
  {
    numeric_index_type lid = static_cast<numeric_index_type>(it);
    auto it_col = find(_gcols.begin(), _gcols.end(), gu2i[it]);
    if (it_col != _gcols.end())
    {
      numeric_index_type gid = std::distance(_gcols.begin(), it_col);
      col_id_cond_u2i_mp.insert(std::make_pair(lid, gid));
    }
    else
      mooseError("DOF ", gu2i[it], " does not exist in the columns of the condensed system");
  }

  for (auto it : index_range(gu2c))
  {
    numeric_index_type lid = static_cast<numeric_index_type>(it);
    auto it_col = find(_gcols.begin(), _gcols.end(), gu2c[it]);
    if (it_col != _gcols.end())
    {
      numeric_index_type gid = std::distance(_gcols.begin(), it_col);
      col_id_cond_u2c_mp.insert(std::make_pair(lid, gid));
    }
    else
      mooseError("DOF ", gu2c[it], " does not exist in the columns of the condensed system");
  }

  // if ((!row_id_cond_mp.empty()) && (!col_id_cond_u2i_mp.empty()))
  // {
  _J_condensed->add_sparse_matrix(*MDinvK2ci, row_id_cond_u2i_mp, col_id_cond_u2i_mp, -1.0);
  _J_condensed->close();
  // }
#ifdef DEBUG
  std::cout << "row_id_cond_u2i_mp = ";
  for (auto i : row_id_cond_u2i_mp)
    std::cout << "(" << i.first << "," << i.second << ") ";
  std::cout << std::endl;

  std::cout << "col_id_cond_u2i_mp = ";
  for (auto i : col_id_cond_u2i_mp)
    std::cout << "(" << i.first << "," << i.second << ") ";
  std::cout << std::endl;
#endif

#ifdef DEBUG
  std::cout << "Norms of _J_condensed after adding MDinvK2cc: l1-norm = " << _J_condensed->l1_norm()
            << "; infinity-norm = " << _J_condensed->linfty_norm()<<"\n";
            // << "; Frobenius-norm = " << _J_condensed->frobenius_norm() << std::endl;

  // _J_condensed->print_personal();
#endif

  // if ((!row_id_cond_mp.empty()) && (!col_id_cond_u2c_mp.empty()))
  // {
  _J_condensed->add_sparse_matrix(*MDinvK2cc, row_id_cond_u2c_mp, col_id_cond_u2c_mp, -1.0);
  _J_condensed->close();
  // }

#ifdef DEBUG
  std::cout << "row_id_cond_u2c_mp = ";
  for (auto i : row_id_cond_u2c_mp)
    std::cout << "(" << i.first << "," << i.second << ") ";
  std::cout << std::endl;

  std::cout << "col_id_cond_u2c_mp = ";
  for (auto i : col_id_cond_u2c_mp)
    std::cout << "(" << i.first << "," << i.second << ") ";
  std::cout << std::endl;
#endif

#ifdef DEBUG
  std::cout << "Norms of _J_condensed after adding MDinvK2cc: l1-norm = " << _J_condensed->l1_norm()
            << "; infinity-norm = " << _J_condensed->linfty_norm()<<"\n";
            // << "; Frobenius-norm = " << _J_condensed->frobenius_norm() << std::endl;

  _J_condensed->print_personal();
  // _matrix->print_matlab("J_original.mat");
  // _J_condensed->print_matlab("J_sorted.mat");
#endif
}

void
ParallelDualMortarPreconditioner::init()
{
  TIME_SECTION(_init_timer);
  if (!_save_dofs)
  {
    // Get DOFs on the secondary/primary and in subdomains for each variable
    getDofVarInterface();
    getDofVarInterior();

    // Get local DOFs on the secondary/primary and in subdomains for each variable
    getLocalDofVarInterface();
    getLocalDofVarInterior();

    // get local row and col dofs for the condensed Jacobian
    std::vector<dof_id_type> u1c = _local_dof_sets_primary[0];
    std::vector<dof_id_type> u2c = _local_dof_sets_secondary[0];

    std::vector<dof_id_type> lm = _local_dof_sets_secondary[1];

    std::vector<dof_id_type> u1i = _local_dof_sets_interior[0][_primary_subdomain];
    std::vector<dof_id_type> u2i = _local_dof_sets_interior[0][_secondary_subdomain];
    // get local row dofs
    _rows.reserve(_dofmap->n_local_dofs() - u2c.size());
    _rows.insert(_rows.end(), u1i.begin(), u1i.end());
    _rows.insert(_rows.end(), u1c.begin(), u1c.end());
    _rows.insert(_rows.end(), u2i.begin(), u2i.end());
    _rows.insert(_rows.end(), lm.begin(), lm.end());
    // sort rows
    std::sort(_rows.begin(), _rows.end());
    // get local col dofs
    _cols.reserve(_dofmap->n_local_dofs() - lm.size());
    _cols.insert(_cols.end(), u1i.begin(), u1i.end());
    _cols.insert(_cols.end(), u1c.begin(), u1c.end());
    _cols.insert(_cols.end(), u2i.begin(), u2i.end());
    // sort cols
    // have gcols in a certain order so that the resultant matrix has no zero-diagonal terms
    std::sort(_cols.begin(), _cols.end());
    _cols.insert(_cols.end(), u2c.begin(), u2c.end());

    // get global row and col dofs for the condensed Jacobian
    u1c = _dof_sets_primary[0];
    u2c = _dof_sets_secondary[0];

    lm = _dof_sets_secondary[1];

    u1i = _dof_sets_interior[0][_primary_subdomain];
    u2i = _dof_sets_interior[0][_secondary_subdomain];

    // get global row dofs
    _grows.reserve(_dofmap->n_dofs() - u2c.size());
    _grows.insert(_grows.end(), u1i.begin(), u1i.end());
    _grows.insert(_grows.end(), u1c.begin(), u1c.end());
    _grows.insert(_grows.end(), u2i.begin(), u2i.end());
    std::sort(_grows.begin(), _grows.end());
    _grows.insert(_grows.end(), lm.begin(), lm.end());

    // get global col dofs
    _gcols.reserve(_dofmap->n_dofs() - lm.size());
    _gcols.insert(_gcols.end(), u1i.begin(), u1i.end());
    _gcols.insert(_gcols.end(), u1c.begin(), u1c.end());
    _gcols.insert(_gcols.end(), u2i.begin(), u2i.end());
    // sort first
    // have gcols in a certain order so that the resultant matrix has no zero-diagonal terms
    std::sort(_gcols.begin(), _gcols.end());
    _gcols.insert(_gcols.end(), u2c.begin(), u2c.end());

    // sort the rows and cols
    _grows_unsrt.resize(_grows.size());
    _grows_unsrt = _grows;
    _gcols_unsrt.resize(_gcols.size());
    _gcols_unsrt = _gcols;

    // vector for mappping old system indices to the new system
    _grow_hat.resize(_dofmap->n_dofs());
    std::fill(_grow_hat.begin(), _grow_hat.end(), libMesh::invalid_uint);
    for (auto i : index_range(_grows))
      _grow_hat[_grows[i]] = i;

    _gcol_hat.resize(_dofmap->n_dofs());
    std::fill(_gcol_hat.begin(), _gcol_hat.end(), libMesh::invalid_uint);
    for (auto i : index_range(_gcols))
      _gcol_hat[_gcols[i]] = i;

    for (auto i : index_range(_dof_sets_secondary[0]))
      _u2c_to_idx.insert(std::make_pair(_dof_sets_secondary[0][i], i));
    // sort u2c indices, need to happen AFTER all grows and gcols have been saved
    std::sort(_local_dof_sets_secondary[0].begin(), _local_dof_sets_secondary[0].end());
    std::sort(_dof_sets_secondary[0].begin(), _dof_sets_secondary[0].end());

#ifdef DEBUG
    std::cout << "_grows = ";
    for (auto i : _grows)
      std::cout << i << " ";
    std::cout << std::endl;
    // std::cout << "_grows_unsrt = ";
    // for (auto i : _grows_unsrt)
    //   std::cout << i << " ";
    // std::cout << std::endl;

    std::cout << "_gcols = ";
    for (auto i : _gcols)
      std::cout << i << " ";
    std::cout << std::endl;
    // std::cout << "_gcols_unsrt = ";
    // for (auto i : _gcols_unsrt)
    //   std::cout << i << " ";
    // std::cout << std::endl;
#endif

    _save_dofs = true;
  }

  if (!_preconditioner)
    _preconditioner =
        Preconditioner<Number>::build_preconditioner(MoosePreconditioner::_communicator);

  _is_initialized = true;
}

void
ParallelDualMortarPreconditioner::print_node_info()
{
  NodeRange * range = _mesh->getActiveNodeRange();
  for (const auto & node : *range)
  {
    node->print_info();
  }
}

void
ParallelDualMortarPreconditioner::setup()
{
  condenseSystem();

  // make sure diagonal entries are not empty
  for (auto i = _J_condensed->row_start(); i < _J_condensed->row_stop(); ++i)
    _J_condensed->add(i, i, 0.0);
  _J_condensed->close();

  _preconditioner->set_matrix(*_J_condensed);
  _preconditioner->set_type(_pre_type);
  _preconditioner->init();
}

void
ParallelDualMortarPreconditioner::apply(const NumericVector<Number> & y, NumericVector<Number> & x)
{
  TIME_SECTION(_apply_timer);

  std::vector<dof_id_type> lm = _dof_sets_secondary[1];

  getCondensedXY(y, x);

#ifdef DEBUG
  std::cout << "Before Apply: \n";
  std::cout << "\tx norm = " << x.l2_norm() << "\n";
  std::cout << "\ty norm = " << y.l2_norm() << "\n";
  std::cout << "\t_x_hat norm = " << _x_hat->l2_norm() << "\n";
  std::cout << "\t_y_hat norm = " << _y_hat->l2_norm() << "\n";
  // std::cout<<"x_hat = \n"<<(*_x_hat)<<std::endl;
  // std::cout<<"y_hat = \n"<<(*_y_hat)<<std::endl;
#endif

  _preconditioner->apply(*_y_hat, *_x_hat);

#ifdef DEBUG
  std::cout << "After Apply: \n";
  std::cout << "\t_x_hat norm  = " << _x_hat->l2_norm() << "\n";
  std::cout << "\t_y_hat norm  = " << _y_hat->l2_norm() << "\n";
  // std::cout<<"x_hat = \n"<<(*_x_hat)<<std::endl;
  // std::cout<<"y_hat = \n"<<(*_y_hat)<<std::endl;
#endif

  computeLM();

#ifdef DEBUG
  std::cout << "\t_lambda norm  = " << _lambda->l2_norm() << "\n";
#endif

  // create a local copy of the vector
  std::unique_ptr<NumericVector<Number>> x_hat_localized(
      NumericVector<Number>::build(MoosePreconditioner::_communicator));
  x_hat_localized->init(_J_condensed->n(), false, SERIAL);

  _x_hat->localize(*x_hat_localized);
  x_hat_localized->close();

  std::unique_ptr<NumericVector<Number>> lambda_localized(
      NumericVector<Number>::build(MoosePreconditioner::_communicator));
  lambda_localized->init(_D->m(), false, SERIAL);

  _lambda->localize(*lambda_localized);
  lambda_localized->close();

  // update x
  for (dof_id_type id1 = 0; id1 < _gcols.size(); ++id1)
  {
    dof_id_type id0 = _gcols[id1]; // id in the original system
    // if (x.is_local(id0))
    if (id0 >= x.first_local_index() && id0 < x.last_local_index())
      x.set(id0, (*x_hat_localized)(id1));
  }

  for (dof_id_type id1 = 0; id1 < lm.size(); ++id1)
  {
    dof_id_type id0 = lm[id1]; // id in the original system
    // if (x.is_local(id0))
    if (id0 >= x.first_local_index() && id0 < x.last_local_index())
      x.set(id0, (*lambda_localized)(id1));
  }

  x.close();

#ifdef DEBUG
  // std::cout << "\tx  = " << x << "\n";
  // std::cout << "\ty  = " << y << "\n";
  std::cout << "\tx_norm  = " << x.l2_norm() << "\n";
  std::cout << "\ty_norm  = " << y.l2_norm() << "\n";
#endif
}

void
ParallelDualMortarPreconditioner::getCondensedXY(const NumericVector<Number> & y,
                                                 NumericVector<Number> & x)
{
  std::vector<dof_id_type> u1c = _dof_sets_primary[0];
  std::vector<dof_id_type> u2c = _dof_sets_secondary[0];

  std::vector<dof_id_type> u2c_local = _local_dof_sets_secondary[0];

  _x_hat = x.zero_clone();
  _y_hat = y.zero_clone();
  _x_hat->init(_J_condensed->n(), _J_condensed->local_n(), false, PARALLEL);
  _y_hat->init(_J_condensed->m(), _J_condensed->local_m(), false, PARALLEL);

  x.create_subvector(*_x_hat, _gcols);
  y.create_subvector(*_y_hat, _grows);

  _r2c = y.zero_clone();
  _r2c->init(_MDinv->n(), _MDinv->local_n(), false, PARALLEL);

  std::unique_ptr<NumericVector<Number>> mdinv_r2c(
      NumericVector<Number>::build(MoosePreconditioner::_communicator));
  mdinv_r2c->init(_MDinv->m(), _MDinv->local_m(), false, PARALLEL);

  std::unique_ptr<NumericVector<Number>> mdinv_r2c_localized(
      NumericVector<Number>::build(MoosePreconditioner::_communicator));
  mdinv_r2c_localized->init(mdinv_r2c->size(), false, SERIAL);

  // get _r2c from the original y
  y.create_subvector(*_r2c, u2c);
  _r2c->close();

  _MDinv->vector_mult(*mdinv_r2c, *_r2c);
  mdinv_r2c->close();

  mdinv_r2c->localize(*mdinv_r2c_localized);
  mdinv_r2c_localized->close();

  for (auto idx : index_range(_grows))
  {
    dof_id_type id0 = _grows[idx]; // row id in the original system
    // if id0 is in u1c, then need to subtract
    // otherwise, copy from y
    auto it_row = find(u1c.begin(), u1c.end(), id0);

    if (it_row != u1c.end())
    {
      // if (_y_hat->is_local(idx))
      if (idx >= _y_hat->first_local_index() && idx < _y_hat->last_local_index())
      {
        Number temp = (*_y_hat)(idx);
        _y_hat->set(idx, temp - (*mdinv_r2c_localized)(std::distance(u1c.begin(), it_row)));
      }
    }
  }

  _y_hat->close();
  _x_hat->close();
}

void
ParallelDualMortarPreconditioner::computeLM()
{
  std::vector<dof_id_type> lm = _local_dof_sets_secondary[1];

  std::vector<dof_id_type> u2i = _dof_sets_interior[0][_secondary_subdomain];
  std::vector<dof_id_type> u2c = _dof_sets_secondary[0];

  std::vector<dof_id_type> u1c = _dof_sets_primary[0];
  std::vector<dof_id_type> u1i = _dof_sets_interior[0][_primary_subdomain];

  _lambda = _r2c->zero_clone();
  _lambda->init(_D->m(), _D->local_m(), false, PARALLEL);

  _x2i = _r2c->zero_clone();
  _x2i->init(_K2ci->n(), _K2ci->local_n(), false, PARALLEL);

  _x2c = _r2c->zero_clone();
  _x2c->init(_K2cc->n(), _K2cc->local_n(), false, PARALLEL);

  // get x2i, x2c from _x_hat
  std::vector<numeric_index_type> x2i_indices, x2c_indices;
  for (auto i : u2i)
  {
    if (_gcol_hat[i] == libMesh::invalid_uint)
      mooseError("Check indices for u2i.");
    x2i_indices.push_back(_gcol_hat[i]);
  }

  for (auto i : u2c)
  {
    if (_gcol_hat[i] == libMesh::invalid_uint)
      mooseError("Check indices for u2c.");
    x2c_indices.push_back(_gcol_hat[i]);
  }

#ifdef DEBUG
  std::cout << "x2i_indices = ";
  for (auto i : x2i_indices)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "x2c_indices = ";
  for (auto i : x2c_indices)
    std::cout << i << " ";
  std::cout << std::endl;
#endif

  _x_hat->create_subvector(*_x2i, x2i_indices);
  _x_hat->create_subvector(*_x2c, x2c_indices);

  _x2i->close();
  _x2c->close();

#ifdef DEBUG
  std::cout << "_x2i norm  = " << _x2i->l2_norm() << std::endl;
  std::cout << "_x2c norm  = " << _x2c->l2_norm() << std::endl;
#endif

  std::unique_ptr<NumericVector<Number>> vec = _r2c->zero_clone(), tmp = _r2c->clone();
  // vec=_K2ci*_x2i;
  _K2ci->vector_mult(*vec, *_x2i);
  (*tmp) -= (*vec);
  vec->close();
  tmp->close();
#ifdef DEBUG
  std::cout << "tmp norm  = " << tmp->l2_norm() << std::endl;
#endif
  // vec=_K2cc*_x2c;
  _K2cc->vector_mult(*vec, *_x2c);
  (*tmp) -= (*vec);
  vec->close();
  tmp->close();
#ifdef DEBUG
  std::cout << "tmp norm  = " << tmp->l2_norm() << std::endl;
#endif
  _D->vector_mult(*_lambda, *tmp);
  _lambda->close();
}

void
ParallelDualMortarPreconditioner::clear()
{
}
