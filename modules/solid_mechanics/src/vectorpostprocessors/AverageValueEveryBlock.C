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
#include "SystemBase.h"
#include "Assembly.h"
#include "libmesh/quadrature.h"

registerMooseObject("SolidMechanicsApp", AverageValueEveryBlock);

InputParameters
AverageValueEveryBlock::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();

  params.addClassDescription("Compute the section's variable average in three-dimensions "
                             "given a user-defined definition of the cross section.");
  // params.addParam<std::vector<SubdomainName>>(
  //     "block",
  //     "The list of blocks in which to search for cross sectional nodes to compute the variable "
  //     "average.");
  // params.addRequiredParam<Point>("axis_direction", "Direction of the structural component's
  // axis"); params.addRequiredParam<Point>("reference_point",
  //                                "Structural component reference starting point from which the "
  //                                "input parameter 'lengths' applies.");
  // params.addParam<Real>("cross_section_maximum_radius",
  //                       std::numeric_limits<double>::max(),
  //                       "Radial distance with respect to the body axis within which nodes are "
  //                       "considered to belong to this "
  //                       "structural component. Used to disambiguate multiple components that
  //                       share " "the same mesh block.");

  params.addRequiredParam<std::vector<VariableName>>("variables",
                                                     "Variables for block-averaged output.");

  // params.addRequiredParam<std::vector<Real>>(
  //     "lengths",
  //     "Distance(s) along axis of from reference_point at which to compute average values.");
  // params.addParam<Real>("tolerance",
  //                       1.0e-6,
  //                       "Maximum axial distance of nodes from the specified axial lengths to "
  //                       "consider them in the cross-section average");
  return params;
}

AverageValueEveryBlock::AverageValueEveryBlock(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    // _mesh(_fe_problem.mesh()),
    _mesh(_subproblem.mesh()),
    _assembly(_subproblem.assembly(0, _sys.number())),
    _q_point(_assembly.qPoints()),
    _qrule(_assembly.qRule()),
    _JxW(_assembly.JxW()),
    _coord(_assembly.coordTransformation()),
    _variables(getParam<std::vector<VariableName>>("variables")),
    _num_cols(_variables.size() + 1), // add one colum for the subdomain ID
    _num_rows(_mesh.meshSubdomains().size())
// _direction(getParam<Point>("axis_direction")),
// _reference_point(getParam<Point>("reference_point")),
// _lengths(getParam<std::vector<Real>>("lengths")),
// _tolerance(getParam<Real>("tolerance")),
// _number_of_nodes(_lengths.size()),
// _cross_section_maximum_radius(getParam<Real>("cross_section_maximum_radius"))
{
  _output_vector.resize(_num_cols);
  for (const auto j : make_range(_num_cols))
  {
    if (j == 0)
      _output_vector[j] = &declareVector("subdomain_id");
    else
      _output_vector[j] = &declareVector(_variables[j - 1]);
  }

  // for (const auto j : make_range(_variables.size()))
  // {
  //   const MooseVariable & variable = _sys.getFieldVariable<Real>(_tid, _variables[j]);
  //   if (!variable.isNodal())
  //     paramError(
  //         "variables",
  //         "The variables provided to this vector postprocessor must be defined at the nodes.");
  // }
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

  _block_ids = _mesh.meshSubdomains();
  int idx = 0;
  for (const auto sid : _block_ids)
    _block_id_map[sid] = idx++;
}

void
AverageValueEveryBlock::finalize()
{
  for (const auto row : make_range(_num_rows))
    for (const auto col : make_range(_num_cols))
    {
      if (!MooseUtils::absoluteFuzzyEqual(_subdomain_areas[col], 0.0))
        (*_output_vector[col])[row] /= _subdomain_areas[col];
    }
}

void
AverageValueEveryBlock::execute()
{
  // Get the current solution vectors
  const NumericVector<Number> & sol = *_sys.currentSolution();

  for (const auto & elem : *_mesh.getActiveLocalElementRange())
  {
    auto sid = elem->subdomain_id();

    // Prepare the element for integration
    _fe_problem.setCurrentSubdomainID(elem, 0);
    _fe_problem.prepare(elem, 0);
    _fe_problem.reinitElem(elem, 0);

    // Calculate the integral of specified variables inside the element
    for (const auto j : make_range(_variables.size()))
    {
      const MooseVariable & var = _sys.getFieldVariable<Real>(_tid, _variables[j]);
      // auto integral = computeIntegral(var.sln());

      Real sum = 0;

      for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
        sum += _JxW[qp] * _coord[qp] * var.sln()[qp];

      // return sum;
    }

    auto idx = _block_id_map[sid];

    // std::cout << "subdomain ID: " << sid << "; " << "element ID: " << elem_ptr->id() <<
    // std::endl;

    // calculate the integral of variables inside this element

    std::cout << "current element ID: " << _assembly.elem()->id() << std::endl;
  }
}

// Real
// AverageValueEveryBlock::computeIntegral(const FieldVariableValue & sol) const
// {
//   Real sum = 0;

//   for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
//     sum += _JxW[qp] * _coord[qp] * sol[qp];

//   return sum;
// }

// Real
// AverageValueEveryBlock::distancePointToPlane(const Node & node,
//                                              const Point & reference_point,
//                                              const Real length) const
// {
//   // Compute node location w.r.t. structural component length
//   const Point relative_distance{
//       node(0) - reference_point(0), node(1) - reference_point(1), node(2) - reference_point(2)};

//   const Real axial_distance = _direction * relative_distance;
//   const Real in_plane_distance =
//       (relative_distance - relative_distance * _direction * _direction).norm();

//   // If the in-plane distance is greater than the specified cross-section radius, the point is not
//   // in this component
//   if (in_plane_distance > _cross_section_maximum_radius)
//     return std::numeric_limits<double>::max();

//   return std::abs(axial_distance - length);
// }
