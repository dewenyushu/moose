# DMP

The Dual Mortar Preconditioner (DMP) is designed to allow the usage of scalable iterative solvers/preconditioners for the system of equations that have saddle point characteristics. The saddle point type of Jacobian often comes from the enforcement of constraints using Lagrange multipliers. Its special numerical character prevents the usage of many scalable iterative solvers. To resolve this issue, a static condensation step is carried out in DMP to remove the degree of freedoms that are associated with the Lagrange multipliers. This results in a positive definite system which can be solved by a broader range of solvers/preconditioners with improved efficiency. Note, the DMP is computationally advantageous when the +dual basis+ is utilized. One can enable the usage of +dual basis+ by enabling `use_dual = true` in the `Variables` block:

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

## Example Input File Syntax

!listing moose/test/tests/preconditioners/dmp/dmp_test.i block=Preconditioning


This initial implementation of DMP works for `EqualValueConstraint` and requires the information about the primary and secondary subdomains. Making it compatible for more complicated problems, such as mortar-based thermal and mechanical contact, remains an ongoing effort.

!syntax parameters /Preconditioning/DMP

!syntax inputs /Preconditioning/DMP

!syntax children /Preconditioning/DMP
