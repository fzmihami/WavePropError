#pragma once

/**
 * @file Dispersion.h
 * @brief Header file defining the Dispersion class for computing dispersive
 * terms in the numerical wave solver.
 *
 * The Dispersion class is responsible for calculating the dispersive terms in
 * the continuity and momentum equations.
 */

#include "GlobalVariables.h"
#include "TimeVariables.h"
#include <iostream>

/**
 * @class Dispersion
 * @brief Stores and computes dispersive terms for the numerical wave solver.
 *
 * The `Dispersion` class calculates and applies dispersive terms in the
 * continuity and momentum equations. It also provides functionality for solving
 * a tridiagonal system using the Thomas algorithm.
 */
class Dispersion {
private:
  unsigned int N; ///< Number of grid points in the computational domain

public:
  float *PhiC; ///< Dispersive term for the continuity equation
  float *PhiM; ///< Dispersive term for the momentum equation
  float *Uxx;  ///< Second-order derivative of the velocity field

  float *D1; ///< Lower diagonal terms of the tridiagonal matrix
  float *D2; ///< Upper diagonal terms of the tridiagonal matrix
  float *D3; ///< Diagonal terms of the tridiagonal matrix

public:
  /**
   * @brief Constructor for the Dispersion class.
   *
   * Initializes the dispersion terms and precomputes the diagonal terms of the
   * tridiagonal matrix.
   *
   * @param S0 Reference to a State object containing wave field variables at
   * the current time step.
   */
  Dispersion(const State &S0);

  /**
   * @brief Computes the dispersive terms for the continuity and momentum
   * equations.
   *
   * This function calculates the second-order derivative of the velocity field
   * and updates the dispersive terms.
   *
   * @param S0 Reference to a State object containing wave field variables at
   * the current time step.
   */
  void ComputeDispersion(State &S0);

  /**
   * @brief Solves the tridiagonal matrix using the Thomas algorithm.
   *
   * This function computes the velocity variable by solving the tridiagonal
   * system.
   *
   * @param S1 Reference to a State object containing wave field variables at
   * the current time step.
   */
  void SolveTridiagonalMatrix(State &S1);

  /**
   * @brief Destructor for the Dispersion class.
   *
   * Cleans up dynamically allocated memory for the dispersion terms and
   * tridiagonal matrix.
   */
  ~Dispersion();
};

/**
 * @brief Solves a tridiagonal system of equations using the Thomas algorithm.
 *
 * This function modifies the input arrays `d` in place to contain the solution
 * to the system.
 *
 * @param a Lower diagonal coefficients of the tridiagonal matrix
 * @param b Main diagonal coefficients of the tridiagonal matrix
 * @param c Upper diagonal coefficients of the tridiagonal matrix
 * @param d Right-hand side vector
 * @param n Size of the system
 */
void ThomasAlgorithm(const float *a, const float *b, const float *c, float *d,
                     int n);