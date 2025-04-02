/**
 * @file Boundaries.cpp
 * @brief Implementation of functions for applying the sponge-layer boundary
 * condition in the numerical wave solver to absorb incoming waves.
 *
 * This file contains the implementation of the `ApplySpongeLayer` function,
 * which applies a sponge layer to the left and right boundaries of the
 * computational domain. The sponge layer is used to absorb incoming waves and
 * prevent reflections from the boundaries.
 */

#include <algorithm>
#include <cstring>
#include <iostream>
#include <math.h>
#include <omp.h>

#include "Boundaries.h"

void ApplySpongeLayer(State &S1) {
  unsigned int N = S1.N;

// left sponge layer
#pragma omp parallel for
  for (unsigned int i = 0; i < nSP; i++) {
    float Cs = pow(alphaSP, pow(gammaSP, i));

    S1.Eta[i] = S1.Eta[i] / Cs;
    S1.Hn[i] = S1.Eta[i] - S1.Zn[i];

    S1.Pn[i] = S1.Pn[i] / Cs;
    S1.Un[i] = S1.Un[i] / Cs;
  }

  // right sponge layer
  if (S1.Grid == grid::staggered) {
#pragma omp parallel for
    for (unsigned int i = N - nSP; i < N; i++) {
      float Cs = pow(alphaSP, pow(gammaSP, N - 1 - i));

      S1.Eta[i] = S1.Eta[i] / Cs;
      S1.Hn[i] = S1.Eta[i] - S1.Zn[i];

      S1.Pn[i + 1] = S1.Pn[i + 1] / Cs;
      S1.Un[i + 1] = S1.Un[i + 1] / Cs;
    }
  } else {
#pragma omp parallel for
    for (unsigned int i = N - nSP; i < N; i++) {
      float Cs = pow(alphaSP, pow(gammaSP, N - 1 - i));

      S1.Eta[i] = S1.Eta[i] / Cs;
      S1.Hn[i] = S1.Eta[i] - S1.Zn[i];

      S1.Pn[i] = S1.Pn[i] / Cs;
      S1.Un[i] = S1.Un[i] / Cs;
    }
  }
}