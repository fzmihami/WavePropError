#pragma once

/**
 * @file SolveEquation.h
 * @brief Implements functions for solving the conservative form of Nwogu's equations.
 *
 * This file provides the implementation of solver routines for Nwogu's equations,
 * supporting both staggered and collocated grid configurations. It includes logic
 * for handling various time integration methods and numerical flux schemes.
 */

#include "Dispersion.h"
#include "FluxComputation.h"
#include "ReadUserInput.h"
#include "TimeVariables.h"
#include "Wavemaker.h"

/**
 * @brief Solves the Nwogu's equation using a staggered grid for the first step
 * of the Runge-Kutta time integration.
 *
 * This function computes the fluxes and source terms for the continuity and
 * momentum equations, and updates the state variables accordingly. It is the
 * first step in the Runge-Kutta time integration.
 *
 * @param S0 The state variable at the current time step.
 * @param S1 The state variable at the next time step.
 * @param In The input parameters for the simulation.
 * @param Fstaggered The flux computation object for the staggered grid.
 * @param Disp The dispersion computation object.
 * @param WM The wavemaker object for source terms.
 * @param dt1 The time step size for the simulation.
 */
void SolveEquationStaggered_step1(State &S0, State &S1, const Input &In,
                                  flux_staggered &Fstaggered, Dispersion &Disp,
                                  WaveMaker &WM, const float &dt1);

/**
 * @brief Solves the Nwogu's equation using a staggered grid for the second step
 * of the Runge-Kutta time integration.
 *
 * This function computes the fluxes and source terms for the continuity and
 * momentum equations, and updates the state variables accordingly. It is the
 * second step in the Runge-Kutta time integration.
 *
 * @param S0 The state variable at the current time step.
 * @param S1 The state variable at the next time step.
 * @param S2 The state variable at the next time step after RK2.
 * @param In The input parameters for the simulation.
 * @param Fstaggered The flux computation object for the staggered grid.
 * @param Disp The dispersion computation object.
 * @param WM The wavemaker object for source terms.
 * @param dt1 The time step size for the simulation.
 */
void SolveEquationStaggered_step2(State &S0, State &S1, State &S2,
                                  const Input &In, flux_staggered &Fstaggered,
                                  Dispersion &Disp, WaveMaker &WM,
                                  const float &dt1);

/**
 * @brief Solves the Nwogu's equation using a staggered grid for the third step
 * of the Runge-Kutta time integration.
 *
 * This function computes the fluxes and source terms for the continuity and
 * momentum equations, and updates the state variables accordingly. It is the
 * third step in the Runge-Kutta time integration.
 *
 * @param S0 The state variable at the current time step.
 * @param S2 The state variable at the next time step after RK2.
 * @param S3 The state variable at the next time step after RK3.
 * @param In The input parameters for the simulation.
 * @param Fstaggered The flux computation object for the staggered grid.
 * @param Disp The dispersion computation object.
 * @param WM The wavemaker object for source terms.
 * @param dt1 The time step size for the simulation.
 */
void SolveEquationStaggered_step3(State &S0, State &S2, State &S3,
                                  const Input &In, flux_staggered &Fstaggered,
                                  Dispersion &Disp, WaveMaker &WM,
                                  const float &dt1);

/**
 * @brief Solves the Nwogu's equation using a collocated grid for the first step
 * of the Runge-Kutta time integration.
 *
 * This function computes the fluxes and source terms for the continuity and
 * momentum equations, and updates the state variables accordingly. It is the
 * first step in the Runge-Kutta time integration.
 *
 * @param S0 The state variable at the current time step.
 * @param S1 The state variable at the next time step.
 * @param In The input parameters for the simulation.
 * @param Fcollocated The flux computation object for the collocated grid.
 * @param Disp The dispersion computation object.
 * @param WM The wavemaker object for source terms.
 * @param dt1 The time step size for the simulation.
 */
void SolveEquationCollocated_step1(State &S0, State &S1, const Input &In,
                                   flux &Fcollocated, Dispersion &Disp,
                                   WaveMaker &WM, const float &dt1);

/**
 * @brief Solves the Nwogu's equation using a collocated grid for the second
 * step of the Runge-Kutta time integration.
 *
 * This function computes the fluxes and source terms for the continuity and
 * momentum equations, and updates the state variables accordingly. It is the
 * second step in the Runge-Kutta time integration.
 *
 * @param S0 The state variable at the current time step.
 * @param S1 The state variable at the next time step.
 * @param S2 The state variable at the next time step after RK2.
 * @param In The input parameters for the simulation.
 * @param Fcollocated The flux computation object for the collocated grid.
 * @param Disp The dispersion computation object.
 * @param WM The wavemaker object for source terms.
 * @param dt1 The time step size for the simulation.
 */
void SolveEquationCollocated_step2(State &S0, State &S1, State &S2,
                                   const Input &In, flux &Fcollocated,
                                   Dispersion &Disp, WaveMaker &WM,
                                   const float &dt1);

/**
 * @brief Solves the Nwogu's equation using a collocated grid for the third step
 * of the Runge-Kutta time integration.
 *
 * This function computes the fluxes and source terms for the continuity and
 * momentum equations, and updates the state variables accordingly. It is the
 * third step in the Runge-Kutta time integration.
 *
 * @param S0 The state variable at the current time step.
 * @param S2 The state variable at the next time step after RK2.
 * @param S3 The state variable at the next time step after RK3.
 * @param In The input parameters for the simulation.
 * @param Fcollocated The flux computation object for the collocated grid.
 * @param Disp The dispersion computation object.
 * @param WM The wavemaker object for source terms.
 * @param dt1 The time step size for the simulation.
 */
void SolveEquationCollocated_step3(State &S0, State &S2, State &S3,
                                   const Input &In, flux &Fcollocated,
                                   Dispersion &Disp, WaveMaker &WM,
                                   const float &dt1);
