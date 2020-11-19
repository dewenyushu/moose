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
  params.addRequiredParam<std::string>("variable",
                                       "Name of the variable that is to be condensed out. Usually "
                                       "this will be the Lagrange multiplier variable.");
  return params;
}

ParallelDualMortarPreconditioner::ParallelDualMortarPreconditioner(const InputParameters & params)
  : MoosePreconditioner(params),
    Preconditioner<Number>(MoosePreconditioner::_communicator),
    _nl(_fe_problem.getNonlinearSystemBase()),
    _mesh(&_fe_problem.mesh()),
    _dofmap(&_nl.system().get_dof_map()),
    _n_vars(_nl.nVariables()),
    _var_name(getParam<std::string>("variable")),
    _var_id(_nl.system().has_variable(_var_name) ? _nl.system().variable_number(_var_name)
                                                 : libMesh::invalid_uint),
    _D(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _M(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _MDinv(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _u2c_rows(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _J_condensed(libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator)),
    _x_hat(NumericVector<Number>::build(MoosePreconditioner::_communicator)),
    _y_hat(NumericVector<Number>::build(MoosePreconditioner::_communicator)),
    _r2c(NumericVector<Number>::build(MoosePreconditioner::_communicator)),
    _lambda(NumericVector<Number>::build(MoosePreconditioner::_communicator)),
    _primary_boundary(getParam<BoundaryID>("primary_boundary")),
    _secondary_boundary(getParam<BoundaryID>("secondary_boundary")),
    _primary_subdomain(getParam<SubdomainID>("primary_subdomain")),
    _secondary_subdomain(getParam<SubdomainID>("secondary_subdomain")),
    _save_dofs(false),
    _init_timer(registerTimedSection("init", 2)),
    _apply_timer(registerTimedSection("apply", 1))
{
  // check if SubdomainID & BoundaryID & variable name are valid
  if (_mesh->meshSubdomains().find(_secondary_subdomain) == _mesh->meshSubdomains().end())
    mooseError("secondary subdomain ID ", _secondary_subdomain, " does not exist.");
  if (_mesh->meshSubdomains().find(_primary_subdomain) == _mesh->meshSubdomains().end())
    mooseError("primary subdomain ID ", _primary_subdomain, " does not exist.");
  if (_mesh->getBoundaryIDs().find(_secondary_boundary) == _mesh->getBoundaryIDs().end())
    mooseError("Secondary boundary ID ", _secondary_boundary, " does not exist.");
  if (_mesh->getBoundaryIDs().find(_primary_boundary) == _mesh->getBoundaryIDs().end())
    mooseError("Secondary boundary ID ", _primary_boundary, " does not exist.");
  if (_var_id == libMesh::invalid_uint)
    paramError("variable", "variable does not exist in the system.");

#ifdef DEBUG
  std::cout << "LM variable ID = " << _var_id << std::endl;
  std::cout << "LM variable name = " << _var_name << std::endl;
#endif

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
ParallelDualMortarPreconditioner::getDofToCondense()
{
  NodeRange * active_nodes = _mesh->getActiveNodeRange();
  for (const auto & node : *active_nodes)
  {
    std::vector<dof_id_type> di;
    _dofmap->dof_indices(node, di, _var_id);
    for (auto index : di)
    {
      _glm.push_back(index);
      if (_dofmap->local_index(index))
        _lm.push_back(index);
    }
  }
  std::sort(_glm.begin(), _glm.end());
  std::sort(_lm.begin(), _lm.end());
}

void
ParallelDualMortarPreconditioner::getDofContact()
{
  for (unsigned int vn = 0; vn < _n_vars; vn++)
  {
    // exclude the lagrange multiplier dofs
    if (vn == _var_id)
      continue;
    // loop over boundary nodes
    ConstBndNodeRange & range = *_mesh->getBoundaryNodeRange();
    std::vector<dof_id_type> di;
    for (const auto & bnode : range)
    {
      const Node * node_bdry = bnode->_node;
      BoundaryID boundary_id = bnode->_bnd_id;

      if (boundary_id == _primary_boundary)
      {
        _dofmap->dof_indices(node_bdry, di, vn);
        for (auto index : di)
        {
          _gu1c.push_back(index);
          if (_dofmap->local_index(index))
            _u1c.push_back(index);
        }
      }

      if (boundary_id == _secondary_boundary)
      {
        _dofmap->dof_indices(node_bdry, di, vn);
        // get corresponding dof of lm on the secondary boundary
        std::vector<dof_id_type> di_lm;
        _dofmap->dof_indices(node_bdry, di_lm, _var_id);
        libmesh_assert(di_lm.size() == di.size());
        for (auto i : index_range(di))
        {
          _gu2c.push_back(di[i]);
          _map_glm_gu2c.insert(std::make_pair(di_lm[i], di[i]));
          _map_gu2c_glm.insert(std::make_pair(di[i], di_lm[i]));
          if (_dofmap->local_index(di[i]))
            _u2c.push_back(di[i]);
        }
      }
    }
  }

  std::sort(_gu1c.begin(), _gu1c.end());
  std::sort(_u1c.begin(), _u1c.end());
  std::sort(_gu2c.begin(), _gu2c.end());
  std::sort(_u2c.begin(), _u2c.end());

  for (auto i : index_range(_glm))
    _map_gu2c_to_order.insert(std::make_pair(_map_glm_gu2c.at(_glm[i]), i));
}

void
ParallelDualMortarPreconditioner::getDofColRow()
{
  // row: all without u2c
  // col: all without lm
  for (numeric_index_type i = 0; i < _dofmap->n_dofs(); ++i)
  {
    if (std::find(_gu2c.begin(), _gu2c.end(), i) != _gu2c.end())
      continue;

    _grows.push_back(i);
    if (_dofmap->local_index(i))
      _rows.push_back(i);

    // ensure the lm and u2c correspondance
    if (_map_glm_gu2c.find(i) != _map_glm_gu2c.end())
    {
      auto u2c_idx = _map_glm_gu2c.at(i);
      _gcols.push_back(u2c_idx);

      if (_dofmap->local_index(u2c_idx))
        _cols.push_back(u2c_idx);
    }
    else
    {
      _gcols.push_back(i);

      if (_dofmap->local_index(i))
        _cols.push_back(i);
    }
  }
}

void
ParallelDualMortarPreconditioner::init()
{
  TIME_SECTION(_init_timer);
  if (!_save_dofs)
  {
    // save necessary dofs
    getDofToCondense();
    getDofContact();

    // get condensed dofs for rows and cols
    getDofColRow();

#ifdef DEBUG
    std::cout << "_glm = ";
    for (auto i : _glm)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_lm = ";
    for (auto i : _lm)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_gu1c = ";
    for (auto i : _gu1c)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_u1c = ";
    for (auto i : _u1c)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_gu2c = ";
    for (auto i : _gu2c)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_u2c = ";
    for (auto i : _u2c)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_map_glm_gu2c = ";
    for (auto i : _map_glm_gu2c)
      std::cout << "(" << i.first << ", " << i.second << ") \n";
    std::cout << std::endl;

    std::cout << "_grows = ";
    for (auto i : _grows)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_rows = ";
    for (auto i : _rows)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_gcols = ";
    for (auto i : _gcols)
      std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "_cols = ";
    for (auto i : _cols)
      std::cout << i << " ";
    std::cout << std::endl;


    std::cout << "_map_gu2c_to_order = ";
    for (auto i : _map_gu2c_to_order)
      std::cout << "(" << i.first << ", " << i.second << ") \n";
    std::cout << std::endl;
#endif

    _save_dofs = true;
  }

  if (!_preconditioner)
    _preconditioner =
        Preconditioner<Number>::build_preconditioner(MoosePreconditioner::_communicator);

  _is_initialized = true;
}

void
ParallelDualMortarPreconditioner::condenseSystem()
{
#ifdef DEBUG
  std::cout << "_u1c = ";
  for (auto i : _u1c)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "_lm = ";
  for (auto i : _lm)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "_glm = ";
  for (auto i : _glm)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "_u2c = ";
  for (auto i : _u2c)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "_gu2c = ";
  for (auto i : _gu2c)
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

  _matrix->create_submatrix(*_M, _u1c, _lm);

#ifdef DEBUG
  // std::cout << "_M = \n";
  // _M->print_personal();
  std::cout << "Norms of _M_transpose: l1-norm = " << _M->l1_norm()
            << "; infinity-norm = " << _M->linfty_norm() << "\n";
  // << "; frobenius_norm = " << _M->frobenius_norm() << std::endl;
#endif

  // get the row associated with u2c
  _u2c_rows->init(_gu2c.size(), _gcols.size(), _u2c.size(), _cols.size());
  MatSetOption(_u2c_rows->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  _matrix->create_submatrix_nosort(*_u2c_rows, _gu2c, _gcols);

#ifdef DEBUG
  std::cout << "Norms of _u2c_rows: l1-norm = " << _u2c_rows->l1_norm()
            << "; infinity-norm = " << _u2c_rows->linfty_norm() << "\n";
#endif

  // invert _D:
  // _D should be strictly diagonal if dual_mortar approach is utilized
  // so we only need to compute the reciprocal number of the diagonal entries
  // to save memory, no new matrix is created
  _matrix->create_submatrix(*_D, _u2c, _lm);
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
    diag_D->set(_map_gu2c_to_order.at(_gu2c[i]), (*_D)(i, _map_gu2c_to_order.at(_gu2c[i])));

  _D->zero();

  for (numeric_index_type i = _D->row_start(); i < _D->row_stop(); ++i)
    if (!MooseUtils::absoluteFuzzyEqual((*diag_D)(i), 0.0))
      _D->set(i, _map_gu2c_to_order.at(_gu2c[i]), 1.0 / (*diag_D)(i));
  _D->close();

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
            << _J_condensed->l1_norm() << "; infinity-norm = " << _J_condensed->linfty_norm()
            << "\n";
  // << "; frobenius_norm = " << _J_condensed->frobenius_norm() << std::endl;
  // _J_condensed->print_personal();
#endif

  // compute changed parts: MDinv* [K2ci, K2cc]
  std::unique_ptr<PetscMatrix<Number>> MDinv_u2c_rows(
      libmesh_make_unique<PetscMatrix<Number>>(MoosePreconditioner::_communicator));

  _MDinv->matrix_matrix_mult(*_u2c_rows, *MDinv_u2c_rows);

#ifdef DEBUG
  std::cout << "Norms of MDinv_u2c_rows: l1-norm = " << MDinv_u2c_rows->l1_norm()
            << "; infinity-norm = " << MDinv_u2c_rows->linfty_norm() << "\n";
#endif

  // add changed parts to _J_condensed
  // original system row_id: u1c
  // original system col_id: _gcols
  // std::vector<numeric_index_type> row_id_cond, col_id_cond_u2i, col_id_cond_u2c;
  std::map<numeric_index_type, numeric_index_type> row_id_mp, col_id_mp;

  for (auto it : index_range(_gu1c))
  {
    numeric_index_type lid = static_cast<numeric_index_type>(it);
    auto it_row = find(_grows.begin(), _grows.end(), _gu1c[it]);
    if (it_row != _grows.end())
    {
      numeric_index_type gid = std::distance(_grows.begin(), it_row);
      if (lid >= MDinv_u2c_rows->row_start() && lid < MDinv_u2c_rows->row_stop())
        row_id_mp.insert(std::make_pair(lid, gid));
    }
    else
      mooseError("DOF ", _gu1c[it], " does not exist in the rows of the condensed system");
  }

  // map for cols
  for (auto it : index_range(_gcols))
  {
    numeric_index_type lid = static_cast<numeric_index_type>(it);
    col_id_mp.insert(std::make_pair(lid, lid));
  }

  // one-step add-matrix
  _J_condensed->add_sparse_matrix(*MDinv_u2c_rows, row_id_mp, col_id_mp, -1.0);
  _J_condensed->close();

#ifdef DEBUG
  std::cout << "Norms of _J_condensed after adding MDinvK2cc: l1-norm = " << _J_condensed->l1_norm()
            << "; infinity-norm = " << _J_condensed->linfty_norm() << "\n";
  // << "; Frobenius-norm = " << _J_condensed->frobenius_norm() << std::endl;

  // _J_condensed->print_personal();
#endif
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

  getFullSolution(y, x);

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
  // create a copy of y
  std::unique_ptr<NumericVector<Number>> y_copy(
      NumericVector<Number>::build(MoosePreconditioner::_communicator));
  y_copy = y.clone();

  _x_hat->init(_J_condensed->n(), _J_condensed->local_n(), false, PARALLEL);
  _y_hat->init(_J_condensed->m(), _J_condensed->local_m(), false, PARALLEL);

  x.create_subvector(*_x_hat, _gcols);

  _r2c->init(_MDinv->n(), _MDinv->local_n(), false, PARALLEL);

  std::unique_ptr<NumericVector<Number>> mdinv_r2c(
      NumericVector<Number>::build(MoosePreconditioner::_communicator));
  mdinv_r2c->init(_MDinv->m(), _MDinv->local_m(), false, PARALLEL);

  // get _r2c from the original y
  y.create_subvector(*_r2c, _gu2c);

  _MDinv->vector_mult(*mdinv_r2c, *_r2c);
  mdinv_r2c->close();

  // change values in y_copy
  std::vector<numeric_index_type> dof_indices;
  std::vector<Number> vals;
  for (auto idx = mdinv_r2c->first_local_index(); idx < mdinv_r2c->last_local_index(); ++idx)
  {
    dof_indices.push_back(_gu1c[idx]);
    vals.push_back(-(*mdinv_r2c)(idx)); // note the minus sign here
  }

  y_copy->add_vector(vals.data(), dof_indices);

  y_copy->create_subvector(*_y_hat, _grows);

  _y_hat->close();
  _x_hat->close();
}

void
ParallelDualMortarPreconditioner::computeLM()
{
  _lambda->init(_D->m(), _D->local_m(), false, PARALLEL);

  std::unique_ptr<NumericVector<Number>> u2c_rows_x_hat(
      NumericVector<Number>::build(MoosePreconditioner::_communicator));
  u2c_rows_x_hat->init(_u2c_rows->m(), _u2c_rows->local_m(), false, PARALLEL);
  _u2c_rows->vector_mult(*u2c_rows_x_hat, *_x_hat);
  u2c_rows_x_hat->close();

  (*_r2c) -= (*u2c_rows_x_hat);
  _r2c->close();
  _D->vector_mult(*_lambda, *_r2c);
  _lambda->close();
}

void
ParallelDualMortarPreconditioner::getFullSolution(const NumericVector<Number> & /*y*/,
                                                  NumericVector<Number> & x)
{
  std::vector<numeric_index_type> dof_indices;
  std::vector<Number> vals;

  // save values and indices from _x_hat and _lambda
  for (auto i = _x_hat->first_local_index(); i < _x_hat->last_local_index(); ++i)
  {
    dof_indices.push_back(_gcols[i]);
    vals.push_back((*_x_hat)(i));
  }

  for (auto i = _lambda->first_local_index(); i < _lambda->last_local_index(); ++i)
  {
    dof_indices.push_back(_glm[i]);
    vals.push_back((*_lambda)(i));
  }

  x.insert(vals.data(), dof_indices);
  x.close();
}

void
ParallelDualMortarPreconditioner::clear()
{
}
