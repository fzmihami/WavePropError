/**
 * @file GlobalVariables.cpp
 * @brief This file contains the global constants and variables used in the
 * numerical solver for wave propagation.
 *
 * The global variables are used throughout the model to define physical
 * constants, numerical scheme settings, and sponge layer parameters.
 *
 * Constants
 * ---------
 * G : float
 *     Gravitational acceleration [m/s^2].
 * PI : float
 *     Mathematical constant pi.
 * hmin : float
 *     Minimum water depth threshold used in computations to avoid small values
 * of the total water depth. alpha1 : float Constant used in Nwogu's equations
 * to achieve optimized dispersion properties for kh in the range [0, pi].
 *     default: -0.39
 * gamma1 : float
 *     Constant used in Nwogu's equations to improve dispersion accuracy.
 *     default: -0.531
 * alphaSP : float
 *     Damping factor used in the sponge layer calculation for wave energy
 * dissipation. default: 2 gammaSP : float Scaling factor for sponge layer's
 * damping effect, typically in the range [0.88, 0.92]. default: 0.92 nSP :
 * unsigned int Number of grid points in the sponge layer on each boundary.
 *     default: 100 (Range: [50, 100]).
 * epsilonSP : float
 *     Threshold value above which the sponge layer is considered inactive.
 *     default: 0.999
 * thetaRC : float
 *     Parameter for Reconstruction of Shallow Water Equations (SWE).
 *     default: 1.0
 * positionWM : unsigned int
 *     Position of the wave maker within the model grid.
 *     default: 0
 *     The correct position of the wave maker is set in the Wavemaker.cpp file.
 *
 * Variables
 * ---------
 * ListeScheme : std::string[4]
 *     List of numerical schemes available for simulation:
 *     - "conservative_staggered": Conservative staggered grid scheme.
 *     - "hllc": Harten-Lax-van Leer Contact (HLLC) scheme.
 *     - "central_upwind": Central upwind scheme.
 *     - "hll": Harten-Lax-van Leer (HLL) scheme.
 *
 * @note The global variables are defined in this file and are used throughout
 * the numerical solver (C++ codes).
 */

#include "GlobalVariables.h"

const float G = 9.81;
const float PI = 3.14159265358979323846;
float hmin = 1e-6;

const std::string ListeScheme[4] = {"conservative_staggered", "hllc",
                                    "central_upwind", "hll"};

const float alpha1 = -0.39; // Nwogu
const float gamma1 = -0.531;

float alphaSP = 2;      // value suggested 2
float gammaSP = 0.92;   // Range [0.88 0.92]
unsigned int nSP = 100; // Range [50 - 100]
float epsilonSP =
    0.999; // for spl>epsilonSP we consider that the sponge layer is not active

const float thetaRC = 1.0; // Reconstruction SWE

unsigned positionWM = 0; // position of the wavemaker
