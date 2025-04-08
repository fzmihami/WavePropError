#pragma once

/**
 * @file GlobalVariables.h
 * @brief Defines global constants and variables used in the wave propagation solver.
 *
 * This header provides physical constants, numerical parameters, and sponge-layer
 * settings shared across the numerical wave model.
 *
 * ## Constants
 * - `G` : Gravitational acceleration \f$ [m/s^2] \f$
 * - `PI` : Mathematical constant pi
 * - `hmin` : Minimum water depth to prevent dry cells
 * - `alpha1`, `gamma1` : Coefficients from Nwogu's equations for dispersion optimization
 * - `thetaRC` : Reconstruction parameter in the shallow water scheme
 *
 * ## Sponge Layer Parameters
 * - `alphaSP` : Damping factor
 * - `gammaSP` : Damping strength scaling factor, typically in [0.88, 0.92]
 * - `nSP` : Number of grid points for sponge layer per boundary
 * - `epsilonSP` : Threshold above which the sponge layer becomes inactive
 *
 * ## Numerical Schemes
 * - `ListeScheme` : List of available numerical flux schemes
 *   - "conservative_staggered"
 *   - "hllc"
 *   - "central_upwind"
 *   - "hll"
 *
 * ## Wave Generation
 * - `positionWM` : Index of the wavemaker location; typically set in `Wavemaker.cpp`
 *
 * @note The global variables are defined in this file and are used throughout
 * the numerical solver (C++ codes).
 */

 #include <iostream>

 // Grid and scheme enumerations
 enum scheme { conservative_staggered, hllc, central_upwind, hll };
 enum grid { collocated, staggered };
 
 // List of numerical schemes
 extern const std::string ListeScheme[4];
 
 // Physical constants
 extern const float G;
 extern const float PI;
 extern float hmin;
 
 // Nwogu’s dispersion coefficients
 extern const float alpha1;
 extern const float gamma1;
 
 // Sponge layer parameters
 extern float alphaSP;
 extern float gammaSP;
 extern unsigned int nSP;
 extern float epsilonSP;
 
 // Reconstruction parameter for SWE
 extern const float thetaRC;
 
 // Wavemaker position
 extern unsigned int positionWM;