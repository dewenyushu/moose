# VCP

!syntax description /Preconditioning/VCP

## Overview

The Variable Condensation Preconditioner (VCP) is designed to condense out variables from the linear system of equations and apply the preconditioner/solver on the simplified equation. This is useful for problems that have saddle point characteristics, which prohibits the usage of scalable iterative solvers/preconditioners.

 The saddle point type of Jacobian often comes from the enforcement of constraints using Lagrange multipliers. Its special numerical character prevents the usage of many scalable iterative solvers. To resolve this issue, a static condensation step is carried out in VCP to remove the degree of freedoms that are associated with the Lagrange multipliers. This results in a positive definite system which can be solved by a broader range of solvers/preconditioners with improved efficiency. Note, the VCP is computationally advantageous when the +dual basis+ is utilized. One can enable the usage of +dual basis+ by enabling `use_dual = true` in the `Variables` block:

```
[Variables]
  [./lm]
    order = FIRST
    family = LAGRANGE
    block = slave_lower
    use_dual = true
  [../]
[]
```

Or setting the `use_dual = true` in the `Contact` block for mortar-based contact:

```
[Contact]
  [leftright]
    mesh = combined_mesh
    secondary = '11'
    primary = '23'

    use_dual = true

    formulation = mortar
    model = frictionless
  [../]
[]
```

## Example Input File Syntax

!listing moose/test/tests/preconditioners/vcp/vcp_test.i block=Preconditioning


This initial implementation of DMP allows usage of +boomeramg+ for `EqualValueConstraint` and has been tested to work for small-scale mortar-based mechanical and thermal-mechanical contact problems.

!syntax parameters /Preconditioning/VCP

!syntax inputs /Preconditioning/VCP

!syntax children /Preconditioning/VCP
