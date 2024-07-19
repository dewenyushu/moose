//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AverageValueEveryBlock.h"
#include "MooseMesh.h"
#include "Assembly.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "libmesh/quadrature.h"

registerMooseObject("SolidMechanicsApp", AverageValueEveryBlock);

InputParameters
AverageValueEveryBlock::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();

  params.addClassDescription("Compute the section's variable average in three-dimensions "
                             "given a user-defined definition of the cross section.");
  params.addRequiredCoupledVar("variables", "Variables for block-averaged output.");
  return params;
}

AverageValueEveryBlock::AverageValueEveryBlock(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    Coupleable(this, false),
    MooseVariableDependencyInterface(this),
    _mesh(_subproblem.mesh()),
    _assembly(_subproblem.assembly(0, _sys.number())),
    _q_point(_assembly.qPoints()),
    _qrule(_assembly.qRule()),
    _JxW(_assembly.JxW()),
    _coord(_assembly.coordTransformation()),
    _variables(coupledNames("variables")),
    _variable_vals(coupledValues("variables")),
    _num_cols(_variables.size() + 1), // add one colum for the subdomain ID
    _num_rows(_mesh.meshSubdomains().size())
{
  _output_vector.resize(_num_cols);
  for (const auto j : make_range(_num_cols))
  {
    if (j == 0)
      _output_vector[j] = &declareVector("subdomain_id");
    else
      _output_vector[j] = &declareVector(_variables[j - 1]);
  }
}

void
AverageValueEveryBlock::initialize()
{
  for (const auto j : make_range(_num_cols))
  {
    _output_vector[j]->clear();
    _output_vector[j]->resize(_num_rows, 0.0);
  }
  _subdomain_areas.resize(_num_rows, 0.0);

  int idx = 0;
  for (const auto sid : _mesh.meshSubdomains())
  {
    (*_output_vector[0])[idx] = sid;
    _block_id_map[sid] = idx++;
  }
}

void
AverageValueEveryBlock::finalize()
{
  // paralle communication
  for (const auto row : make_range(_num_rows))
    _communicator.sum(_subdomain_areas[row]);

  for (const auto row : make_range(_num_rows))
  {
    for (const auto col : make_range(_num_cols))
    {
      if (col == 0)
        _communicator.max((*_output_vector[col])[row]);
      else
        _communicator.sum((*_output_vector[col])[row]);
    }
  }

  // calculate outputs
  for (const auto row : make_range(_num_rows))
  {
    for (const auto col : make_range(_num_cols))
    {
      if (col == 0) // skip the column for the subdomain ID
        continue;

      if (!MooseUtils::absoluteFuzzyEqual(_subdomain_areas[row], 0.0))
        (*_output_vector[col])[row] /= _subdomain_areas[row];
    }
  }
}

void
AverageValueEveryBlock::execute()
{
  for (const auto & elem : *_mesh.getActiveLocalElementRange())
  {
    // Prepare the element for integration
    _fe_problem.setCurrentSubdomainID(elem, 0);
    _fe_problem.prepare(elem, 0);
    _fe_problem.reinitElem(elem, 0);

    // compute integral of every coupled variable and add to corresponding items in _output_vector
    // compute element area and add to the corresponding entry in _subdomain_areas
    computeIntegral();
  }
}

void
AverageValueEveryBlock::computeIntegral()
{
  auto sid = _assembly.elem()->subdomain_id();

  auto row_id = _block_id_map.at(sid);

  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
  {
    for (const auto i : make_range(_variables.size()))
    {
      (*_output_vector[i + 1])[row_id] += _JxW[qp] * _coord[qp] * (*_variable_vals[i])[qp];
    }

    _subdomain_areas[row_id] += _JxW[qp] * _coord[qp];
  }
}
