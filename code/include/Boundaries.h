#pragma once

/**
 * @file Boundaries.h
 * @brief Header file defining functions for applying the sponge-layer boundary
 * condition in the numerical wave solver to absorb incoming waves.
 */

#include <iostream>

#include "GlobalVariables.h"
#include "TimeVariables.h"

/**
 * @brief Applies the sponge-layer boundary condition to the left and right
 * boundaries of the computational domain. This fucntion is based on the
 * approach described in Larsen and Dancy (1983).
 *
 * This function modifies the following variables in the `State` object `S1` for
 * grid points inside the sponge layer:
 * - `S1.Eta`: The water surface elevation.
 * - `S1.Hn`: The total water depth.
 * - `S1.Un`: The horizontal velocity.
 * - `S1.Pn`: The conserved momentum in Nwogu's equation.
 *
 * @param S1 Reference to the `State` object containing wave field variables at
 * the current time step t.
 */
void ApplySpongeLayer(State &S1);