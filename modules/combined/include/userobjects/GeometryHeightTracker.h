#pragma once

#include <cfloat>

#include "ElementUserObject.h"
#include "Function.h"
#include <cfloat>

class GeometryHeightTracker : public ElementUserObject
{
public:
  static InputParameters validParams();

  GeometryHeightTracker(const InputParameters & parameters);

  Real getLayerHeight() const { return _height; };

protected:
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & /*uo*/) override{};

private:
  void updateHeightForElement(const Elem * elem);

  void update_auxiliary_boundary_info();

  void setHeatSourceCenter()
  {
    _center(0) = _px[0];
    _center(1) = _py[0];
    _center(2) = _pz[0];
  };

  /// center of the heat source center along three directions
  const VariableValue & _px;
  const VariableValue & _py;
  const VariableValue & _pz;

  /// The name of the moving boundary
  const BoundaryName & _moving_boundary_name;

  /// The Id of the moving boundary
  BoundaryID _moving_boundary_id;

  /// The name of the auxiliary boundary
  const BoundaryName & _auxiliary_boundary_name;

  /// The ID of the auxiliary boundary
  BoundaryID _aux_boundary_id;

  /// initial center of the printing path
  std::vector<Real> _initial_center;

  /// previous step center on the printing path (used to check if a new layer is added)
  Point _prev_center;

  /// current center on the printing path
  Point _center;

  /// height for the new center
  Real _height;

  /// flag for starting a new layer
  bool _new_layer;
};
