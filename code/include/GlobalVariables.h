#pragma once

/**
 * @file GlobalVariables.h
 * @brief Header file defining enums, constants, and global variables for the
 * numerical wave solver.
 */

#include <iostream>

enum scheme { conservative_staggered, hllc, central_upwind, hll };
enum grid { collocated, staggered };

extern const std::string ListeScheme[4];

extern const float G;
extern const float PI;
extern float hmin;

// Nwogu's equations
extern const float alpha1;
extern const float gamma1;

// sponge layer
extern float alphaSP;    // value suggested 2
extern float gammaSP;    // Range [0.88 0.92]
extern unsigned int nSP; // Range [50 - 100]
extern float epsilonSP;  // for spl>epsilonSP we consider that the sponge layer
                         // is not active anymore

// Reconstruction SWE
extern const float thetaRC;

extern unsigned positionWM; // position of the wavemaker