//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GmshBoxGenerator.h"

registerMooseObject("MooseApp", GmshBoxGenerator);

InputParameters
GmshBoxGenerator::validParams()
{
  InputParameters params = GmshGeneratorBase::validParams();
  params.set<int>("dim") = 3;
  params.addParam<Real>("x", 1.0, "X direction length");
  params.addParam<Real>("y", 1.0, "Y direction length");
  params.addParam<Real>("z", 1.0, "Z direction length");
  params.addParam<Real>("spacing", 0.2, "Mesh spacing");
  return params;
}

GmshBoxGenerator::GmshBoxGenerator(const InputParameters & parameters)
  : GmshGeneratorBase(parameters),
    _x(getParam<Real>("x")),
    _y(getParam<Real>("y")),
    _z(getParam<Real>("z"))
{
}

void
GmshBoxGenerator::generateGeometry()
{
  Real lc = getParam<Real>("spacing");
  int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
  int p2 = gmsh::model::geo::addPoint(_x, 0, 0, lc);
  int p3 = gmsh::model::geo::addPoint(_x, _y, 0, lc);
  int p4 = gmsh::model::geo::addPoint(0, _y, 0, lc);
  int p5 = gmsh::model::geo::addPoint(0, 0, _z, lc);
  int p6 = gmsh::model::geo::addPoint(_x, 0, _z, lc);
  int p7 = gmsh::model::geo::addPoint(_x, _y, _z, lc);
  int p8 = gmsh::model::geo::addPoint(0, _y, _z, lc);

  int l1 = gmsh::model::geo::addLine(p1, p2);
  int l2 = gmsh::model::geo::addLine(p2, p3);
  int l3 = gmsh::model::geo::addLine(p3, p4);
  int l4 = gmsh::model::geo::addLine(p4, p1);

  int l5 = gmsh::model::geo::addLine(p5, p6);
  int l6 = gmsh::model::geo::addLine(p6, p7);
  int l7 = gmsh::model::geo::addLine(p7, p8);
  int l8 = gmsh::model::geo::addLine(p8, p5);

  int l9 = gmsh::model::geo::addLine(p1, p5);
  int l10 = gmsh::model::geo::addLine(p2, p6);
  int l11 = gmsh::model::geo::addLine(p3, p7);
  int l12 = gmsh::model::geo::addLine(p4, p8);

  int cl1 = gmsh::model::geo::addCurveLoop({l1, l2, l3, l4});
  int cl2 = gmsh::model::geo::addCurveLoop({-l5,  -l8, -l7, -l6});
  int cl3 = gmsh::model::geo::addCurveLoop({-l3, l11, l7, -l12});
  int cl4 = gmsh::model::geo::addCurveLoop({-l1, l9, l5, -l10});
  int cl5 = gmsh::model::geo::addCurveLoop({-l4, l12, l8, -l9});
  int cl6 = gmsh::model::geo::addCurveLoop({-l2, l10, l6, -l11});

  int s1 = gmsh::model::geo::addPlaneSurface({cl1});
  int s2 = gmsh::model::geo::addPlaneSurface({cl2});
  int s3 = gmsh::model::geo::addPlaneSurface({cl3});
  int s4 = gmsh::model::geo::addPlaneSurface({cl4});
  int s5 = gmsh::model::geo::addPlaneSurface({cl5});
  int s6 = gmsh::model::geo::addPlaneSurface({cl6});

  int sl = gmsh::model::geo::addSurfaceLoop({s1, s3, s2, s4, s5, s6});

  gmsh::model::geo::addVolume({sl}, 1);
}
