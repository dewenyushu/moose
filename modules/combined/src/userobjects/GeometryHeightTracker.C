#include "GeometryHeightTracker.h"

registerMooseObject("SMARTAMApp", GeometryHeightTracker);

InputParameters
GeometryHeightTracker::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addRequiredParam<std::vector<Real>>("initial_center",
                                             "The initial position of the heat source center.");

  params.addRequiredCoupledVar(
      "position_x", "Coupled variable representing the x-coordinate of the heat source center.");

  params.addRequiredCoupledVar(
      "position_y", "Coupled variable representing the y-coordinate of the heat source center.");

  params.addRequiredCoupledVar(
      "position_z", "Coupled variable representing the z-coordinate of the heat source center.");

  params.addRequiredParam<BoundaryName>(
      "moving_boundary_name",
      "Boundary to modify when an element is moved. A boundary with the provided name will be "
      "created if not already exists on the mesh.");

  params.addParam<BoundaryName>(
      "auxiliary_boundary_name", "auxiliary_boundary", "Boundary from the previous layer.");

  params.addClassDescription(
      "Element subdomain modifier that tracks the height of the previous layer and set the heat "
      "source center on the previous layer interface.");

  return params;
}

GeometryHeightTracker::GeometryHeightTracker(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _px(coupledValue("position_x")),
    _py(coupledValue("position_y")),
    _pz(coupledValue("position_z")),
    _moving_boundary_name(getParam<BoundaryName>("moving_boundary_name")),
    _auxiliary_boundary_name(getParam<BoundaryName>("auxiliary_boundary_name")),
    _initial_center(getParam<std::vector<Real>>("initial_center")),
    _new_layer(false)
{
  // get info of the tracking boundary
  _moving_boundary_id = _mesh.getBoundaryIDs({{_moving_boundary_name}}, true)[0];

  if (_moving_boundary_id == Moose::INVALID_BOUNDARY_ID)
    paramError("moving_boundary_name", "Invalid moving boundary name.");

  // Get the BoundaryID for the auxiliary boundary from the mesh
  _aux_boundary_id = _mesh.getBoundaryIDs({{_auxiliary_boundary_name}}, true)[0];

  if (_initial_center.size() != 3)
    paramError("initial_center",
               "Three components are required for setting the initial center of the heat source.");

  _prev_center = Point(_initial_center[0], _initial_center[1], _initial_center[2]);
  _center = _prev_center;
}

void
GeometryHeightTracker::initialize()
{
  // update the heat source center position
  setHeatSourceCenter();

  // if a new layer starts printing, we need to reset the element list
  if (_center(2) - _prev_center(2) > TOLERANCE)
    _new_layer = true;

  // we would like to find the maximum height, so set it to minimum first
  _height = -DBL_MAX;
}

void
GeometryHeightTracker::execute()
{
}

void
GeometryHeightTracker::update_auxiliary_boundary_info()
{
  _mesh.getMesh().get_boundary_info().remove_id(_aux_boundary_id);

  // Get the BoundaryID for the auxiliary boundary and set proper name
  _aux_boundary_id = _mesh.getBoundaryIDs({{_auxiliary_boundary_name}}, true)[0];

  _mesh.getMesh().get_boundary_info().sideset_name(_aux_boundary_id) = _auxiliary_boundary_name;
  _mesh.getMesh().get_boundary_info().nodeset_name(_aux_boundary_id) = _auxiliary_boundary_name;

  // loop through elements in the boundary
  auto & elem_side_bnd_ids = _mesh.getMesh().get_boundary_info().get_sideset_map();

  // add sides in moving boundary to the aux boundary
  for (const auto & [elem, side_bnd] : elem_side_bnd_ids)
  {
    auto side = side_bnd.first;
    auto boundary_id = side_bnd.second;
    if (boundary_id == _moving_boundary_id)
      _mesh.getMesh().get_boundary_info().add_side(elem, side, _aux_boundary_id);
  }

  _mesh.getMesh().get_boundary_info().parallel_sync_side_ids();
  _mesh.getMesh().get_boundary_info().parallel_sync_node_ids();
  _mesh.update();
  _mesh.getMesh().prepare_for_use();
}

void
GeometryHeightTracker::finalize()
{
  if (_new_layer == true)
  {
    update_auxiliary_boundary_info();
    _console << COLOR_MAGENTA << "Update aux boundary info " << COLOR_DEFAULT << std::endl;
  }

  // loop through elements in the boundary
  auto & elem_side_bnd_ids = _mesh.getMesh().get_boundary_info().get_sideset_map();

  for (const auto & [elem, side_bnd] : elem_side_bnd_ids)
  {
    auto boundary_id = side_bnd.second;
    if (boundary_id == _aux_boundary_id)
      updateHeightForElement(elem);
  }

  // to be MPI safe, gather the maximum value across processors
  gatherMax(_height);

  // for the 1st layer, make sure we do not go to -DBL_MAX in height
  _height = _height < _initial_center[2] ? _initial_center[2] : _height;

  _console << COLOR_MAGENTA << "Adjusted laser center elevation = " << _height << COLOR_DEFAULT
           << std::endl;

  // update the center
  _prev_center = _center;

  // update flag
  _new_layer = false;
}

void
GeometryHeightTracker::updateHeightForElement(const Elem * elem)
{
  // find the bounding coordinate values
  Point max_coord(-DBL_MAX, -DBL_MAX, -DBL_MAX), min_coord(DBL_MAX, DBL_MAX, DBL_MAX);
  for (auto idx : elem->node_index_range())
  {
    auto pt = elem->point(idx);
    for (unsigned int j = 0; j < 3; ++j)
    {
      max_coord(j) = std::max(pt(j), max_coord(j));
      min_coord(j) = std::min(pt(j), min_coord(j));
    }
  }

  Point dummy_center(_center(0), _center(1), (max_coord(2) + min_coord(2)) / 2.0);

  if (elem->contains_point(dummy_center))
    _height = std::max(_height, max_coord(2));
}
