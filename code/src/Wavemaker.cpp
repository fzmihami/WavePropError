/**
 * @file Wavemaker.cpp
 * @brief Implements the WaveMaker class for generating wave forcing through a source term.
 *
 * This file provides the implementation of the WaveMaker class, which computes
 * the source term used to generate waves within the computational domain. It also
 * manages and stores the relevant wavemaker parameters. The wave forcing is
 * constructed by superimposing contributions from multiple frequency components.
 */

#include <algorithm>
#include <cstring>
#include <iostream>
#include <math.h>
#include <omp.h>

#include "Wavemaker.h"

WaveMaker::WaveMaker(const Input &In, const unsigned int &N) {
  this->Nf = In.WaveInput.size();
  this->N = N;
  this->Beta = 2.0f * M_PI / In.LengthDomain;
  this->depth = In.WaterDepth;

  this->Freq = new float[Nf];
  memset(this->Freq, 0, Nf * sizeof(float));
  this->Amp = new float[Nf];
  memset(this->Amp, 0, Nf * sizeof(float));
  this->Phase = new float[Nf];
  memset(this->Phase, 0, Nf * sizeof(float));

  this->Omg = new float[Nf];
  memset(this->Omg, 0, Nf * sizeof(float));
  this->K = new float[Nf];
  memset(this->K, 0, Nf * sizeof(float));
  this->Di = new float[Nf];
  memset(this->Di, 0, Nf * sizeof(float));

  this->SourceTerm = new float[N];
  memset(this->SourceTerm, 0, N * sizeof(float));

  for (unsigned int i = 0; i < Nf; i++) {
    this->Freq[i] = In.WaveInput[i][0];
    this->Amp[i] = In.WaveInput[i][1];
    this->Phase[i] = In.WaveInput[i][2];
  }

  for (unsigned int i = 0; i < Nf; i++) {
    this->Omg[i] = 2.0f * M_PI * this->Freq[i];
    this->K[i] = wavelength(this->Omg[i], In.WaterDepth);
  }

  // compute peak wavelength
  unsigned int max_amp_index =
      std::max_element(this->Amp, this->Amp + Nf) - this->Amp;
  this->Kp = this->K[max_amp_index];
  float Lp = 2.0f * M_PI / this->Kp;
  this->Beta = 80.0f / (de * de * Lp * Lp);

  // compute the number of grid points for the wavemaker
  this->width = 0.5f * Lp;
  this->Ns = (unsigned int)std::round(this->width / In.dx);

  // wavemaker position at 2*Lp from the left boundary
  this->posWM = (unsigned int)std::round(2.0f * Lp / In.dx);
  positionWM = posWM;

  // sponge layer width - Larsen and Dancy 1983
  nSP = (unsigned int)std::round(
      Lp / In.dx); // the sponge layer width is the peak wavelength
  gammaSP = pow(log(1.0 / epsilonSP) / log(alphaSP),
                1.0 / nSP); // extract correspending gamma value

  // verify if the domain is too small for the wavemaker
  if (In.LengthDomain < 5.0f * Lp) {
    std::cerr << "The domain is too small to fit the wavemaker" << std::endl;
    exit(1);
  }

  // Compute Di
  for (unsigned int i = 0; i < Nf; i++) {
    float li = K[i];
    float Bi = sqrt(PI / Beta) * exp(-(li * li / (4 * Beta)));

    Di[i] = (2 * Amp[i] *
             (Omg[i] * Omg[i] - alpha3 * G * pow(K[i], 4) * pow(depth, 3))) /
            (Omg[i] * K[i] * Bi * (1.0 - alpha1 * pow(K[i] * depth, 2)));
  }
}

WaveMaker::~WaveMaker() {
  delete[] Freq;
  delete[] Amp;
  delete[] Phase;
  delete[] Omg;
  delete[] K;
  delete[] Di;
  delete[] SourceTerm;
}

void WaveMaker::ComputeSourceTerm(const float &t, const float &dx) {
  float St = 0.0f;

#pragma omp parallel for reduction(+ : St)
  for (unsigned int i = 0; i < Nf; i++) {
    float ramp = tanh(t * Omg[i] / (2 * PI));

    St += ramp * Di[i] * sin(-Omg[i] * t + Phase[i]);
  }

// compute the source term
#pragma omp parallel for
  for (unsigned int i = 1; i < Ns; i++) {
    SourceTerm[posWM + i] = St * exp(-Beta * pow(i * dx, 2));
    SourceTerm[posWM - i] = St * exp(-Beta * pow(i * dx, 2));
  }
  SourceTerm[posWM] = St;
}

float wavelength(float omega, float depth) // Compute Kn from Omega
{

  float beta = alpha1 + (float)(1.0 / 3.0);
  float K;

  if (abs(beta) >= 1E-8)
    K = SecondorderEquation(G * beta, -(depth * omega * omega * alpha1 + G),
                            depth * omega * omega);
  else
    K = depth * omega * omega / (depth * omega * omega * alpha1 + G);

  float k = sqrt(K) / depth;

  return k;
}

float SecondorderEquation(float a, float b, float c) {
  float delta = b * b - 4 * a * c;
  float x1 = (-b - sqrt(delta)) / (2.0 * a);

  return x1;
}
