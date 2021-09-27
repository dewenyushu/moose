//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GmshTensionSampleGenerator.h"

registerMooseObject("MooseApp", GmshTensionSampleGenerator);

InputParameters
GmshTensionSampleGenerator::validParams()
{
  InputParameters params = GmshGeneratorBase::validParams();
  params.addParam<Real>(
      "coarse_spacing",
      1e-2,
      "The mesh size of the parts that are away from the hole, which can be coarse meshed.");
  params.addParam<Real>(
      "fine_spacing",
      1e-2,
      "The mesh size of the parts that are close to the hole, which are fine meshed.");
  params.addParam<Real>("half_sample_length", 8, "Half of the overal sample length");
  params.addParam<Real>("half_bridge_height", 0.6, "Half of the sample bridge height");
  params.addParam<Real>("half_sample_height", 2.0, "Half of the height of the gripping part");
  params.addParam<Real>("gripping_length", 4.1, "Length of the gripping part");
  params.addParam<Real>("fillet_radius", 0.3, "Fillet radius");

  params.addParam<Real>("right_circle_x", 5.75, "Right circle centroid x coordinate");
  params.addParam<Real>("right_circle_y", 0, "Right circle centroid y coordinate");
  params.addParam<Real>("right_circle_r", 0.75, "Right circle radius");
  params.addParam<Real>("left_circle_x", -5.75, "Left circle centroid x coordinate");
  params.addParam<Real>("left_circle_y", 0, "Left circle centroid y coordinate");
  params.addParam<Real>("left_circle_r", 0.75, "Left circle radius");

  params.addParam<Real>("scale", 1.0, "Scale factor");
  return params;
}

GmshTensionSampleGenerator::GmshTensionSampleGenerator(const InputParameters & parameters)
  : GmshGeneratorBase(parameters),
    _scale(getParam<Real>("scale")),
    _half_sample_length(getParam<Real>("half_sample_length") * _scale),
    _half_bridge_height(getParam<Real>("half_bridge_height") * _scale),
    _half_sample_height(getParam<Real>("half_sample_height") * _scale),
    _gripping_length(getParam<Real>("gripping_length") * _scale),
    _fillet_radius(getParam<Real>("fillet_radius") * _scale),
    _right_circle_x(getParam<Real>("right_circle_x") * _scale),
    _right_circle_y(getParam<Real>("right_circle_y") * _scale),
    _right_circle_r(getParam<Real>("right_circle_r") * _scale),
    _left_circle_x(getParam<Real>("left_circle_x") * _scale),
    _left_circle_y(getParam<Real>("left_circle_y") * _scale),
    _left_circle_r(getParam<Real>("left_circle_r") * _scale),
    _connection_radius(_half_sample_height - _half_bridge_height - _fillet_radius),
    _half_bridge_length(_half_sample_length - _gripping_length - _connection_radius)
{
  if (_connection_radius <= 0)
    mooseError("The radius of the connection part is not positive.");
  if (_half_bridge_length <= 0)
    mooseError("The bridge length is not positive.");
  if (_half_bridge_length + _connection_radius + 2.0 * _right_circle_r >= _half_sample_length)
    paramError("half_sample_length", "Right hole is out of the sample (x-direction)");
  if (_right_circle_y + _right_circle_r >=
          _half_bridge_height + _connection_radius + _fillet_radius ||
      _right_circle_y - _right_circle_r <=
          -_half_bridge_height - _connection_radius - _fillet_radius)
    paramError("right_circle_y", "Right hole is out of the sample (y-direction)");

  if (_half_bridge_length + _connection_radius + 2.0 * _left_circle_r >= _half_sample_length)
    paramError("half_sample_length", "Left hole is out of the sample");
  if (_left_circle_y + _left_circle_r >=
          _half_bridge_height + _connection_radius + _fillet_radius ||
      _left_circle_y - _left_circle_r <= -_half_bridge_height - _connection_radius - _fillet_radius)
    paramError("left_circle_y", "Left hole is out of the sample (y-direction)");
}

void
GmshTensionSampleGenerator::generateGeometry()
{
  Real lc = getParam<Real>("coarse_spacing") * _scale;
  Real lf = getParam<Real>("fine_spacing") * _scale;

  // top right part
  Real p1x = _half_bridge_length;
  Real p1y = _half_bridge_height;

  Real c1x = p1x;
  Real c1y = p1y + _connection_radius;

  Real p2x = c1x + _connection_radius;
  Real p2y = c1y;

  Real p3x = p2x;
  Real p3y = p2y + _fillet_radius;

  Real p4x = _half_sample_length;
  Real p4y = p3y;

  // top right part
  int p1 = gmsh::model::geo::addPoint(p1x, p1y, 0, lf);
  int p2 = gmsh::model::geo::addPoint(p2x, p2y, 0, lf);
  int p3 = gmsh::model::geo::addPoint(p3x, p3y, 0, lf);
  int p4 = gmsh::model::geo::addPoint(p4x, p4y, 0, lc);
  int c1 = gmsh::model::geo::addPoint(c1x, c1y, 0, lf);
  int m1 = addCurveCenter({p1x, p1y}, {p2x, p2y}, {c1x, c1y});

  // bottom right part
  int p5 = gmsh::model::geo::addPoint(p4x, -p4y, 0, lc);
  int p6 = gmsh::model::geo::addPoint(p3x, -p3y, 0, lf);
  int p7 = gmsh::model::geo::addPoint(p2x, -p2y, 0, lf);
  int p8 = gmsh::model::geo::addPoint(p1x, -p1y, 0, lf);
  int c2 = gmsh::model::geo::addPoint(c1x, -c1y, 0, lf);
  int m2 = addCurveCenter({p2x, -p2y}, {p1x, -p1y}, {c1x, -c1y});

  // bottom left part
  int p9 = gmsh::model::geo::addPoint(-p1x, -p1y, 0, lf);
  int p10 = gmsh::model::geo::addPoint(-p2x, -p2y, 0, lf);
  int p11 = gmsh::model::geo::addPoint(-p3x, -p3y, 0, lf);
  int p12 = gmsh::model::geo::addPoint(-p4x, -p4y, 0, lc);
  int c3 = gmsh::model::geo::addPoint(-c1x, -c1y, 0, lf);
  int m3 = addCurveCenter({-p1x, -p1y}, {-p2x, -p2y}, {-c1x, -c1y});

  // top left part
  int p13 = gmsh::model::geo::addPoint(-p4x, p4y, 0, lc);
  int p14 = gmsh::model::geo::addPoint(-p3x, p3y, 0, lf);
  int p15 = gmsh::model::geo::addPoint(-p2x, p2y, 0, lf);
  int p16 = gmsh::model::geo::addPoint(-p1x, p1y, 0, lf);
  int c4 = gmsh::model::geo::addPoint(-c1x, c1y, 0, lf);
  int m4 = addCurveCenter({-p2x, p2y}, {-p1x, p1y}, {-c1x, c1y});

  // center line
  int p17 = gmsh::model::geo::addPoint(0, p1y, 0, lc);
  int p18 = gmsh::model::geo::addPoint(0, -p1y, 0, lc);

  // int center = gmsh::model::geo::addLine(p17, p18);
  // gmsh::model::addPhysicalGroup(1, {center});

  // lines and curves
  std::vector<int> line_tags_loop;
  line_tags_loop.push_back(gmsh::model::geo::addLine(p17, p1));
  line_tags_loop.push_back(gmsh::model::geo::addCircleArc(p1, c1, m1)); // curve
  line_tags_loop.push_back(gmsh::model::geo::addCircleArc(m1, c1, p2)); // curve
  line_tags_loop.push_back(gmsh::model::geo::addLine(p2, p3));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p3, p4));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p4, p5));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p5, p6));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p6, p7));
  line_tags_loop.push_back(gmsh::model::geo::addCircleArc(p7, c2, m2)); // curve
  line_tags_loop.push_back(gmsh::model::geo::addCircleArc(m2, c2, p8)); // curve
  line_tags_loop.push_back(gmsh::model::geo::addLine(p8, p18));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p18, p9));
  line_tags_loop.push_back(gmsh::model::geo::addCircleArc(p9, c3, m3));  // curve
  line_tags_loop.push_back(gmsh::model::geo::addCircleArc(m3, c3, p10)); // curve
  line_tags_loop.push_back(gmsh::model::geo::addLine(p10, p11));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p11, p12));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p12, p13));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p13, p14));
  line_tags_loop.push_back(gmsh::model::geo::addLine(p14, p15));
  line_tags_loop.push_back(gmsh::model::geo::addCircleArc(p15, c4, m4)); // curve
  line_tags_loop.push_back(gmsh::model::geo::addCircleArc(m4, c4, p16)); // curve
  line_tags_loop.push_back(gmsh::model::geo::addLine(p16, p17));

  // container for curve loop tags
  std::vector<int> tags;
  tags.push_back(gmsh::model::geo::addCurveLoop(line_tags_loop));
  // create holes
  tags.push_back(Moose::gmshAddCircleLoop(
      _right_circle_x, _right_circle_y, _right_circle_r, lf, "right_hole"));
  tags.push_back(
      Moose::gmshAddCircleLoop(_left_circle_x, _left_circle_y, _left_circle_r, lf, "left_hole"));

  // create surface
  gmsh::model::geo::addPlaneSurface(tags, 1);

  if (_dim == 3)
  {
    // extrude to 3D mesh
    double h = 2;
    std::vector<std::pair<int, int>> ov;
    gmsh::model::geo::extrude({{2, 1}}, 0, 0, h, ov, {2}, {1}, true);
    gmsh::model::geo::synchronize();

#ifdef DEBUG
    std::cout << "Finish extrusion" << std::endl;
#endif
  }

  // add physical groups for each entity
  std::vector<std::pair<int, int>> entities;
  gmsh::model::getEntities(entities);
  for (auto e : entities)
  {
    // Dimension and tag of the entity:
    int dim = e.first, tag = e.second;

    // Extract tag for the lower dimensional entities
    if (dim != _dim - 1)
      continue;
#ifdef DEBUG
    std::cout << "Boundary tag " << tag << std::endl;
#endif
    std::vector<int> physicalTags;
    gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);
    if (physicalTags.size() == 0)
      gmsh::model::addPhysicalGroup(dim, {tag});
  }
}

int
GmshTensionSampleGenerator::addCurveCenter(const std::vector<Real> & p1_coord,
                                           const std::vector<Real> & p2_coord,
                                           const std::vector<Real> & c_coord)
{
  Real lc = getParam<Real>("fine_spacing") * _scale;

  if (p1_coord.size() != 2 || p2_coord.size() != 2 || c_coord.size() != 2)
    mooseError("Points defining the curve should have dimension of 2");
  std::vector<Real> N_coord, pc_coord;
  N_coord.resize(2);
  pc_coord.resize(2);

  // radius of the curve
  Real r = std::pow((p1_coord[0] - c_coord[0]), 2) + std::pow((p1_coord[1] - c_coord[1]), 2);

  r = std::sqrt(r);

  // calculate middle point of P1 and P2
  for (int i = 0; i < 2; ++i)
    N_coord[i] = (p1_coord[i] + p2_coord[i]) / 2;

  Real den = std::pow((N_coord[0] - c_coord[0]), 2) + std::pow((N_coord[1] - c_coord[1]), 2);
  Real temp_x = r * (N_coord[0] - c_coord[0]) / std::sqrt(den);
  Real temp_y = r * (N_coord[1] - c_coord[1]) / std::sqrt(den);

  if (p1_coord[0] < p2_coord[0])
  {
    if (c_coord[0] + temp_x < p2_coord[0] && c_coord[0] + temp_x > p1_coord[0])
    {
      pc_coord[0] = c_coord[0] + temp_x;
      pc_coord[1] = c_coord[1] + temp_y;
    }
    else
    {
      pc_coord[0] = c_coord[0] - temp_x;
      pc_coord[1] = c_coord[1] - temp_y;
    }
  }
  else
  {
    if (c_coord[0] + temp_x > p2_coord[0] && c_coord[0] + temp_x < p1_coord[0])
    {
      pc_coord[0] = c_coord[0] + temp_x;
      pc_coord[1] = c_coord[1] + temp_y;
    }
    else
    {
      pc_coord[0] = c_coord[0] - temp_x;
      pc_coord[1] = c_coord[1] - temp_y;
    }
  }
  // check correctness
  if (!MooseUtils::absoluteFuzzyEqual(std::pow((pc_coord[0] - c_coord[0]), 2) +
                                          std::pow((pc_coord[1] - c_coord[1]), 2),
                                      r * r,
                                      1e-6))
    mooseError("Center point coordinate calculation wrong");
  return gmsh::model::geo::addPoint(pc_coord[0], pc_coord[1], 0, lc);
}
