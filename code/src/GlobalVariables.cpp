/**
 * @file GlobalVariables.cpp
 * @brief This file contains the global constants and variables used in the
 * numerical solver for wave propagation.
 *
 * The global variables are used throughout the model to define physical
 * constants, numerical scheme settings, and sponge layer parameters.
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
