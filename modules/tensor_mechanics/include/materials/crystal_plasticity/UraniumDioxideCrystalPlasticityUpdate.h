//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CrystalPlasticityStressUpdateBase.h"

/**
 * UraniumDioxideCrystalPlasticityUpdate uses the multiplicative decomposition of the
 * deformation gradient and solves the PK2 stress residual equation at the
 * intermediate configuration to evolve the material state. The internal
 * variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */

class UraniumDioxideCrystalPlasticityUpdate : public CrystalPlasticityStressUpdateBase
{
public:

  static InputParameters validParams();

  UraniumDioxideCrystalPlasticityUpdate(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   */
  virtual void initQpStatefulProperties() override;

  /**
   * Sets the value of the current and previous substep iteration slip system
   * resistance to the old value at the start of the PK2 stress convergence
   * while loop.
   */
  virtual void setInitialConstitutiveVariableValues() override;

  /**
   * Sets the current slip system resistance value to the previous substep value.
   * In cases where only one substep is taken (or when the first) substep is taken,
   * this method just sets the current value to the old slip system resistance
   * value again.
   */
  virtual void setSubstepConstitutiveVariableValues() override;

  /**
   * Stores the current value of the slip system resistance into a separate
   * material property in case substepping is needed.
   */
  virtual void updateSubstepConstitutiveVariableValues() override;

  virtual bool calculateSlipRate() override;

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  /**
   * The forest dislocation density rate is calculated here (following equations 9 - 10 in the
   * notes)
   */
  virtual void calculateStateVariableEvolutionRateComponent() override;

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;

  virtual void calculateSlipResistance() override;

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  /*
   * Determines if the state variables, e.g. defect densities, have converged
   * by comparing the change in the values over the iteration period.
   */
  virtual bool areConstitutiveStateVariablesConverged() override;

  /*
   * Read parameters (activation energy, drag stress, and the hardening matrix) from the auxiliary
   * input file
   */
  void readDislocationSlipResistanceFileParams();

  /// file name of the text file providing the slip resistance information
  const FileName & _parameter_filename;

  ///@{Parameters used in the slip rate update equation
  const Real & _dG0; // effective thermal activation energy
  const Real & _k;   // Boltzmann constant
  const Real & _T;   // absolute temperature (TODO: coupled variable instead?)
  const Real & _p;   // 0 - 1, 1st exponent parameter related to the obstacles resistance profile
  const Real & _q;   // 1 - 2, 2nd exponent parameter related to the obstacles resistance profile
  const unsigned int & _nu; // attempt (or attack) frequency, which is the number of counts that a
                            // dislocation attempts to overcome an obstacle
  const Real & _lambda_eff; // mean distance between obstacles
  ///@}

  ///@{Parameters used in the hardening mechanism
  std::vector<std::vector<Real>> _K; // the hardening matrix (read from file)
  const Real & _tau0;                // lattice friction or Peierls stress
  const Real & _mu;                  // shear modulus
  const Real & _b;                   // magnitude of the burger's vector
  ///@}

  ///@{Parameters used in the dislocation density evolution
  const Real & _k1;         // constant in the dislocation density evolutino equation
  const Real & _chi;        // interaction parameter
  std::vector<Real> _g;     // normalized activation energy (read from file)
  std::vector<Real> _D;     // drag stress (read from file)
  const Real & _gamma_dot0; // reference shear rate
  const Real & _rho0; // initial dislocation density value
  ///@}

  /// State variables
  MaterialProperty<std::vector<Real>> & _forest_dislocation_density;
  const MaterialProperty<std::vector<Real>> & _forest_dislocation_density_old;

  /**
   * Stores the values of the slip system resistance from the previous substep
   * In classes which use dislocation densities, analogous dislocation density
   * substep vectors will be required.
   */
  std::vector<Real> _previous_substep_forest_dislocation_density;

  /**
   * Caches the value of the current slip system resistance immediately prior
   * to the update of the slip system resistance, and is used to calculate the
   * the slip system resistance increment for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison.
   * In classes which use dislocation densities, analogous dislocation density
   * caching vectors will also be required.
   */
  std::vector<Real> _forest_dislocation_density_before_update;

  /// forest dislocation density evolution component (does not need to be stateful property)
  std::vector<Real> _forest_dislocation_density_rate;

  /// Increment of increased resistance for each slip system
  std::vector<Real> _slip_resistance_increment;

  /**
   * Stores the values of the slip system resistance from the previous substep
   * In classes which use dislocation densities, analogous dislocation density
   * substep vectors will be required.
   */
  std::vector<Real> _previous_substep_slip_resistance;

  /**
   * Caches the value of the current slip system resistance immediately prior
   * to the update of the slip system resistance, and is used to calculate the
   * the slip system resistance increment for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison.
   * In classes which use dislocation densities, analogous dislocation density
   * caching vectors will also be required.
   */
  std::vector<Real> _slip_resistance_before_update;
};
