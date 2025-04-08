/**
 * @file FluxComputation.cpp
 * @brief Implementation of flux computation classes for the Shallow Water Equations 
 *        on both collocated and staggered grids.
 *
 * This file defines the `flux` and `flux_staggered` classes used to compute 
 * fluxes in the Shallow Water Equations. The `flux` class supports collocated 
 * grids using the Finite Volume Godunov method, while the `flux_staggered` 
 * class is designed for staggered grids, implementing a conservative staggered 
 * scheme based on Stelling and Duinmeijer (2003).
 */

#include <cstring>
#include <iostream>
#include <math.h>
#include <omp.h>

#include "FluxComputation.h"
#include "GlobalVariables.h"

// Flux staggered class
// ///////////////////////////////////////////////////////////////////////////
flux_staggered::flux_staggered(unsigned int N) {
  this->N = N;
  WS = new float[N];
  memset(WS, 0, N * sizeof(float));
  HUn = new float[N + 1];
  memset(HUn, 0, (N + 1) * sizeof(float));
  HUUn = new float[N];
  memset(HUUn, 0, N * sizeof(float));
}

flux_staggered::~flux_staggered() {
  delete[] WS;
  delete[] HUn;
  delete[] HUUn;
}

void flux_staggered::ComputeFlux(const State &S0,
                                 const unsigned int &OrderReconstruction) {
  if (OrderReconstruction != 1 && OrderReconstruction != 2) {
    std::cout << "Error: Order of reconstruction not implemented for the "
                 "conservative staggered scheme"
              << std::endl;
    exit(1);
  }

// Compute flux continuity
#pragma omp parallel for
  for (unsigned int i = 0; i < (N + 1); i++) {

    if (i == 0) {
      HUn[i] = S0.Un[i] * S0.Hn[i];
    } else if (i == N) {
      HUn[i] = S0.Un[i] * S0.Hn[i - 1];
    }

    else if (OrderReconstruction == 1 || i == 1 || i == N - 1) {
      HUn[i] = S0.Un[i] * upwind(S0.Un[i], S0.Hn[i - 1], S0.Hn[i]);
    }

    else {
      HUn[i] = S0.Un[i] * upwind(S0.Un[i], S0.Hn[i - 2], S0.Hn[i - 1], S0.Hn[i],
                                 S0.Hn[i + 1]);
    }
  }

// Compute flux momentum
#pragma omp parallel for
  for (unsigned int i = 0; i < N; i++) {

    float pr = 0.5 * (HUn[i] + HUn[i + 1]);
    float ur;

    if (i == 0 || i == N - 1 || OrderReconstruction == 1) {
      ur = upwind(pr, S0.Un[i], S0.Un[i + 1]);
    } else {
      ur = upwind(pr, S0.Un[i - 1], S0.Un[i], S0.Un[i + 1], S0.Un[i + 2]);
    }

    HUUn[i] = pr * ur;
  }

// compute wave speed - Stelling and Duinmeijer (2003)
#pragma omp parallel for
  for (unsigned int i = 0; i < N; i++) {

    if (S0.Hn[i] < hmin) {
      WS[i] = 0;
      continue;
    }

    WS[i] = sqrt(G * S0.Hn[i]) + fabs(HUn[i + 1] + HUn[i]) / (2 * S0.Hn[i]);
  }
}

float upwind(float V, float Q1, float Q2) {
  if (V >= 0)
    return Q1;
  else
    return Q2;
}

float upwind(float V, float Q1, float Q2, float Q3, float Q4) {
  if (V >= 0)
    return (Q2 + 0.5 * minmod(thetaRC * (Q3 - Q2), 0.5 * (Q3 - Q1),
                              thetaRC * (Q2 - Q1)));
  else
    return (Q3 - 0.5 * minmod(thetaRC * (Q4 - Q3), 0.5 * (Q4 - Q2),
                              thetaRC * (Q3 - Q2)));
}

//  Flux class Collocated grid
//  ///////////////////////////////////////////////////////////////////////////

flux::flux(unsigned int N) {
  this->N = N;
  WS = new float[N];
  memset(WS, 0, N * sizeof(float));
  Hn = new float[N];
  memset(Hn, 0, N * sizeof(float));
  Qn = new float[N];
  memset(Qn, 0, N * sizeof(float));
  HL = new float[N - 1];
  memset(HL, 0, (N - 1) * sizeof(float));
  HR = new float[N - 1];
  memset(HR, 0, (N - 1) * sizeof(float));
  QL = new float[N - 1];
  memset(QL, 0, (N - 1) * sizeof(float));
  QR = new float[N - 1];
  memset(QR, 0, (N - 1) * sizeof(float));
  FluxH = new float[N - 1];
  memset(FluxH, 0, (N - 1) * sizeof(float));
  FluxU = new float[N - 1];
  memset(FluxU, 0, (N - 1) * sizeof(float));
}

flux::~flux() {
  delete[] WS;
  delete[] Hn;
  delete[] Qn;
  delete[] HL;
  delete[] HR;
  delete[] QL;
  delete[] QR;
  delete[] FluxH;
  delete[] FluxU;
}

void flux::ComputeFlux(const State &S0, const scheme &Scheme,
                       const unsigned int &OrderReconstruction) {
// Compute Qn = Hn*Un
#pragma omp parallel for
  for (unsigned int i = 0; i < N; i++) {
    Hn[i] = S0.Hn[i];
    Qn[i] = S0.Hn[i] * S0.Un[i];
  }

  // Reconstruction
  switch (OrderReconstruction) {
  case 1:
    Recontruction_Order1(Hn, HL, HR, N);
    Recontruction_Order1(Qn, QL, QR, N);
    break;
  case 2:
    Recontruction_Order2(Hn, HL, HR, N);
    Recontruction_Order2(Qn, QL, QR, N);
    break;
  case 3:
    Recontruction_Order3(Hn, HL, HR, N);
    Recontruction_Order3(Qn, QL, QR, N);
    break;
  case 5:
    Recontruction_Order5(Hn, HL, HR, N);
    Recontruction_Order5(Qn, QL, QR, N);
    break;
  default:
    std::cout << "Error: Order of reconstruction not implemented" << std::endl;
    exit(1);
  }

  // Flux computation
  switch (Scheme) {
  case scheme::central_upwind:
    Flux_CentralUpwind(HL, HR, QL, QR, FluxH, FluxU, N);
    break;
  case scheme::hllc:
    Flux_HLLC(HL, HR, QL, QR, FluxH, FluxU, N);
    break;
  case scheme::hll:
    Flux_HLL(HL, HR, QL, QR, FluxH, FluxU, N);
    break;
  default:
    std::cout << "Error: Scheme not implemented" << std::endl;
    exit(1);
  }

// compute wave speed
#pragma omp parallel for
  for (unsigned int i = 0; i < N; i++) {
    WS[i] = fabs(S0.Un[i]) + sqrt(G * Hn[i]);
  }
}

//  Reconstruction functions + limiters
//  //////////////////////////////////////////////////////////

float min3(float a, float b, float c) { return std::min(std::min(a, b), c); }

float max3(float a, float b, float c) { return std::max(std::max(a, b), c); }

float minmod(float a, float b, float c) {
  if (a > 0 && b > 0 && c > 0) {
    return std::min(a, std::min(b, c));
  } else if (a < 0 && b < 0 && c < 0) {
    return std::max(a, std::max(b, c));
  } else {
    return 0;
  }
}

void Recontruction_Order1(float *Hn, float *Hl, float *Hr, unsigned int N) {
#pragma omp parallel for
  for (unsigned int i = 0; i < N - 1; i++) {
    Hl[i] = Hn[i];
    Hr[i] = Hn[i + 1];
  }
}

void Recontruction_Order2(float *Hn, float *Hl, float *Hr, unsigned int N) {
  float *Hx = (float *)malloc(N * sizeof(float));
  memset(Hx, 0, N * sizeof(float));

#pragma omp parallel for
  for (unsigned int i = 1; i < N - 1; i++) {
    Hx[i] = minmod(thetaRC * (Hn[i] - Hn[i - 1]), 0.5 * (Hn[i + 1] - Hn[i - 1]),
                   thetaRC * (Hn[i + 1] - Hn[i]));
  }

#pragma omp parallel for
  for (unsigned int i = 0; i < N - 1; i++) {
    Hl[i] = Hn[i] + 0.5 * Hx[i];
    Hr[i] = Hn[i + 1] - 0.5 * Hx[i + 1];
  }

  free(Hx);
}

void Recontruction_Order3(float *Hn, float *Hl, float *Hr, unsigned int N) {
  float *rL = (float *)malloc(N * sizeof(float));
  float *rR = (float *)malloc(N * sizeof(float));

  float bL, bR;

#pragma omp parallel for
  for (unsigned int i = 1; i < N - 1; i++) {
    rL[i] = (Hn[i + 1] - Hn[i] + 1E-20) / (Hn[i] - Hn[i - 1] + 1E-20);
    rR[i] = (Hn[i] - Hn[i - 1] + 1E-20) / (Hn[i + 1] - Hn[i] + 1E-20);
  }
  rL[0];
  rL[N - 1];
  rR[0];
  rR[N - 1];

#pragma omp parallel for
  for (unsigned int i = 2; i < N - 3; i++) {

    bL = (1.0 + 2 * rL[i]) / 3.0;
    bR = (1.0 + 2 * rR[i + 1]) / 3.0;

    Hl[i] = Hn[i] + 0.5 * std::max(0.0f, min3(2.0f, 2 * rL[i], bL)) *
                        (Hn[i] - Hn[i - 1]);
    Hr[i] = Hn[i + 1] - 0.5 * std::max(0.0f, min3(2.0f, 2 * rR[i + 1], bR)) *
                            (Hn[i + 2] - Hn[i + 1]);
  }
  Hl[0] = Hn[0];
  Hr[0] = Hn[1];
  Hl[1] = Hn[1];
  Hr[1] = Hn[2];

  Hl[N - 2] = Hn[N - 2];
  Hr[N - 2] = Hn[N - 1];
  Hl[N - 3] = Hn[N - 3];
  Hr[N - 3] = Hn[N - 2];

  free(rL);
  free(rR);
}

void Recontruction_Order5(float *Hn, float *Hl, float *Hr, unsigned int N) {
  float *rL = (float *)malloc(N * sizeof(float));
  float *rR = (float *)malloc(N * sizeof(float));

  float bL, bR;

#pragma omp parallel for
  for (int i = 1; i < N - 1; i++) {
    rL[i] = (Hn[i + 1] - Hn[i] + 1E-20) / (Hn[i] - Hn[i - 1] + 1E-20);
    rR[i] = (Hn[i] - Hn[i - 1] + 1E-20) / (Hn[i + 1] - Hn[i] + 1E-20);
  }
  rL[0];
  rL[N - 1];
  rR[0];
  rR[N - 1];

#pragma omp parallel for
  for (unsigned int i = 2; i < N - 3; i++) {

    bL = (-2.0 / rL[i - 1] + 11.0 + 24.0 * rL[i] - 3.0 * rL[i] * rL[i + 1]) /
         30.0;
    bR =
        (-2.0 / rR[i + 2] + 11.0 + 24.0 * rR[i + 1] - 3.0 * rR[i] * rR[i + 1]) /
        30.0;

    Hl[i] = Hn[i] + 0.5 * std::max(0.0f, min3(2.0, 2 * rL[i], bL)) *
                        (Hn[i] - Hn[i - 1]);
    Hr[i] = Hn[i + 1] - 0.5 * std::max(0.0f, min3(2.0, 2 * rR[i + 1], bR)) *
                            (Hn[i + 2] - Hn[i + 1]);
  }
  Hl[0] = Hn[0];
  Hr[0] = Hn[1];
  Hl[1] = Hn[1];
  Hr[1] = Hn[2];

  Hl[N - 2] = Hn[N - 2];
  Hr[N - 2] = Hn[N - 1];
  Hl[N - 3] = Hn[N - 3];
  Hr[N - 3] = Hn[N - 2];

  free(rL);
  free(rR);
}

// Flux compuation
// ///////////////////////////////////////////////////////////////////////////////
void Flux_CentralUpwind(float *Hl, float *Hr, float *Ql, float *Qr,
                        float *FluxH, float *FluxQ, unsigned int N) {
  float hL, hR, uL, uR;
  float FhL, FhR, FuL, FuR;
  float aL, aR;

#pragma omp parallel for
  for (unsigned int i = 0; i < N - 1; i++) {

    hL = Hl[i];
    uL = Ql[i] / hL;

    FhL = Ql[i];
    FuL = Ql[i] * uL + 0.5 * G * Hl[i] * Hl[i];

    // Right state

    hR = Hr[i];
    uR = Qr[i] / hR;

    FhR = Qr[i];
    FuR = Qr[i] * uR + 0.5 * G * Hr[i] * Hr[i];

    // Local Speeds
    aL = min3(0, uR - sqrt(G * hR), uL - sqrt(G * hL));
    aR = max3(0, uR + sqrt(G * hR), uL + sqrt(G * hL));

    // Central Upwind Algorithm
    FluxH[i] = (aR * FhL - aL * FhR) / (aR - aL) +
               (aR * aL * (Hr[i] - Hl[i])) / (aR - aL);
    FluxQ[i] = (aR * FuL - aL * FuR) / (aR - aL) +
               (aR * aL * (Qr[i] - Ql[i])) / (aR - aL);
  }
}

void Flux_HLL(float *Hl, float *Hr, float *Ql, float *Qr, float *FluxH,
              float *FluxQ, unsigned int N) {
  float hL, hR, uL, uR;
  float FhL, FhR, FuL, FuR;
  float qL, qR;
  float h0;
  float sL, sR;

#pragma omp parallel for
  for (unsigned int i = 0; i < N - 1; i++) {

    // Left state
    hL = Hl[i];
    uL = Ql[i] / Hl[i];
    FhL = hL * uL;
    FuL = hL * uL * uL + G * hL * hL / 2;

    // Right state
    hR = Hr[i];
    uR = Qr[i] / Hr[i];
    FhR = hR * uR;
    FuR = hR * uR * uR + G * hR * hR / 2;

    // Compute h*
    h0 = 0.5 * (hL + hR) -
         0.25 * (uR - uL) * (hL + hR) / (sqrt(G * hL) + sqrt(G * hR));

    // Compute QL, QR
    if (h0 > hL) {
      qL = sqrt(0.5 * (h0 + hL) * h0 / (hL * hL));
    } else {
      qL = 1;
    }

    if (h0 > hR) {
      qR = sqrt(0.5 * (h0 + hR) * h0 / (hR * hR));
    } else {
      qR = 1;
    }

    // Left and right speeds
    sL = uL - qL * sqrt(G * hL);
    sR = uR + qR * sqrt(G * hR);

    // HLL algorithm
    if (0 < sL) {
      FluxH[i] = FhL;
      FluxQ[i] = FuL;
    } else if (sR < 0) {
      FluxH[i] = FhR;
      FluxQ[i] = FuR;
    } else {
      FluxH[i] =
          (sR * FhL - sL * FhR + (sL * sR) * (Hr[i] - Hl[i])) / (sR - sL);
      FluxQ[i] =
          (sR * FuL - sL * FuR + (sL * sR) * (Qr[i] - Ql[i])) / (sR - sL);
    }
  }
}

void Flux_HLLC(float *Hl, float *Hr, float *Ql, float *Qr, float *FluxH,
               float *FluxQ, unsigned int N) {
  float hL, hR, uL, uR;
  float FhL, FhR, FuL, FuR;
  float qL, qR;
  float h0, u0;
  float sL, sR, s0; // s0 is the middle wave speed

  float hStarL, hStarR, qStarL, qStarR;

  float kL, kR;

#pragma omp parallel for
  for (unsigned int i = 0; i < N - 1; i++) {
    // Left state
    hL = Hl[i];
    uL = Ql[i] / Hl[i];
    FhL = hL * uL;
    FuL = hL * uL * uL + G * hL * hL / 2;

    // Right state
    hR = Hr[i];
    uR = Qr[i] / Hr[i];
    FhR = hR * uR;
    FuR = hR * uR * uR + G * hR * hR / 2;

    // Compute h*
    h0 = 0.5 * (hL + hR) -
         0.25 * (uR - uL) * (hL + hR) / (sqrt(G * hL) + sqrt(G * hR));
    u0 =
        0.5 * (uL + uR) - (hR - hL) * (sqrt(G * hL) + sqrt(G * hR)) / (hL + hR);

    // Compute QL, QR
    if (h0 > hL) {
      qL = sqrt(0.5 * (h0 + hL) * h0 / (hL * hL));
    } else {
      qL = 1;
    }

    if (h0 > hR) {
      qR = sqrt(0.5 * (h0 + hR) * h0 / (hR * hR));
    } else {
      qR = 1;
    }

    // Left and right speeds
    sL = uL - qL * sqrt(G * hL);
    sR = uR + qR * sqrt(G * hR);

    // Compute the middle wave speed s0
    s0 = (sL * hR * (uR - sR) - sR * hL * (uL - sL)) /
         (hR * (uR - sR) - hL * (uL - sL));

    // Compute U* L/R
    // coefficients U* = K * [1, s0]
    kL = hL * (sL - uL) / (sL - s0);
    kR = hR * (sR - uR) / (sR - s0);

    hStarL = kL;
    hStarR = kR;
    qStarL = kL * s0;
    qStarR = kR * s0;

    // HLLC algorithm
    if (0 <= sL) {
      FluxH[i] = FhL;
      FluxQ[i] = FuL;
    } else if (sR <= 0) {
      FluxH[i] = FhR;
      FluxQ[i] = FuR;
    } else {
      if (0 <= s0) {
        FluxH[i] = FhL + sL * (hStarL - hL);
        FluxQ[i] = FuL + sL * (qStarL - hL * uL);
      } else {
        FluxH[i] = FhR + sR * (hStarR - hR);
        FluxQ[i] = FuR + sR * (qStarR - hR * uR);
      }
    }
  }
}
