/**
 * @file TimeVariables.cpp
 * @brief Implementation of the State class and functions for managing time
 * variables in the numerical wave solver.
 *
 * The file contains the implementation of the State class, which represents the
 * state variables (i.e., water height, velocity, etc.) at the current time
 * step. It includes methods for updating the state and managing memory. The
 * file also contains functions for printing progress, recording time series
 * data at defined gauges, and saving the simulation results to binary files.
 */

#include <cmath>
#include <iostream>
#include <string.h>

#include "TimeVariables.h"

State::State(const Input &In) {
  dx = In.dx;
  dx_1 = 1.0 / dx;
  N = round(In.LengthDomain / dx) + 1;
  depth = In.WaterDepth;
  if (In.Scheme == scheme::conservative_staggered)
    Grid = grid::staggered;
  else
    Grid = grid::collocated;

  // allocate memory
  Eta = new float[N];
  memset(Eta, 0, N * sizeof(float));
  Hn = new float[N];
  memset(Hn, 0, N * sizeof(float));
  Zn = new float[N];
  memset(Zn, 0, N * sizeof(float));

  if (Grid == grid::staggered) {
    Un = new float[N + 1];
    memset(Un, 0, (N + 1) * sizeof(float));
    Pn = new float[N + 1];
    memset(Pn, 0, (N + 1) * sizeof(float));
  } else {
    Un = new float[N];
    memset(Un, 0, N * sizeof(float));
    Pn = new float[N];
    memset(Pn, 0, N * sizeof(float));
  }

  // initial conditions
  for (unsigned int i = 0; i < N; i++) {
    Hn[i] = depth;
    Zn[i] = -depth;
  }
}

void swap_pointers(float **a, float **b) {
  float *tmp = *a;
  *a = *b;
  *b = tmp;
  tmp = nullptr;
}

void State::Update(State &other) {
  swap_pointers(&Eta, &other.Eta);
  swap_pointers(&Hn, &other.Hn);
  swap_pointers(&Un, &other.Un);
  swap_pointers(&Pn, &other.Pn);
}

State::~State() {
  delete[] Eta;
  delete[] Hn;
  delete[] Un;
  delete[] Pn;
  delete[] Zn;
}

// //////////////////////////////////////////////////////////////////////////////////////////

//// printout to binary File ////////////////////////////
void printbin(float *ptr, int size, std::string namefile1) {
  FILE *fptr;

  const char *namefile = namefile1.c_str();

  if ((fptr = fopen(namefile, "wb")) == NULL) {
    printf("Error! opening file");
    exit(1);
  }
  fwrite(ptr, sizeof(float), size, fptr);
  fclose(fptr);
}

//// Printout Info //////////////////////////////
void printINFO(const Input &In) {
  float info[10];

  info[0] = In.Time;
  info[1] = In.dtOutput;
  info[2] = round(In.Time / In.dtOutput) + 1;

  info[3] = In.dtGauges;
  info[4] = round(In.Time / In.dtGauges) + 1;

  info[5] = In.LengthDomain;
  info[6] = In.dx;
  info[7] = round(In.LengthDomain / In.dx) + 1;

  info[8] = In.WaterDepth;

  info[9] = positionWM;

  std::string namefile = In.ResFolder + "/info.bin";
  printbin(info, 10, namefile);
}

void PrintVariables(State &S0, State &S2, const Input &In, const float &t1,
                    const float &dt1, unsigned int &Kprint, float &Tprint) {
  const unsigned int N = S0.N;

  while (t1 <= Tprint && t1 + dt1 >= Tprint && Tprint <= In.Time) {

    linear_interpolation(S0.Eta, S2.Eta, N, t1, dt1,
                         Tprint); // The intrpolated value is saved in Eta -
                                  // wont be used later Memory recycling
    std::string namefile = In.ResFolder + "/FreeSurface/H_" +
                           intToStringWithLeadingZeros(Kprint, 6) + ".bin";
    printbin(S0.Eta, N, namefile);

    Kprint++;
    Tprint = Kprint * In.dtOutput; // Better in terms of floating precision
  }
}

void RecordTimeSeries(State &S0, State &S2, const Input &In,
                      std::vector<float> *Hgauges, const float &t1,
                      const float &dt1, unsigned int &Kgauges, float &Tgauges) {
  while (t1 <= Tgauges && t1 + dt1 >= Tgauges &&
         (Tgauges <= In.Time || std::abs(Tgauges - In.Time) < 5e-5)) {

    for (unsigned int k = 0; k < In.IndexGauges.size(); k++) {
      int i = In.IndexGauges[k];
      float value =
          linear_interpolation(S0.Eta[i], S2.Eta[i], t1, dt1, Tgauges);
      Hgauges[k].push_back(value);
    }

    Kgauges++;
    Tgauges = Kgauges * In.dtGauges; // Better in terms of floating precision
  }
}

void PrintTimeSeries(std::vector<float> *Hgauges, const Input &In) {
  for (unsigned int k = 0; k < In.IndexGauges.size(); k++) {
    std::string namefile =
        In.ResFolder + "/TimeSeries/Gauge_" + std::to_string(k) + ".bin";
    printbin(&Hgauges[k][0], Hgauges[k].size(), namefile);
  }
}

void linear_interpolation(float *H1, float *H2, unsigned int N, const float &t,
                          const float &dt, const float &Tprint) {
#pragma omp parallel for
  for (unsigned int k = 0; k < N; k++) {
    H1[k] = (H2[k] - H1[k]) * (Tprint - t) / dt + H1[k];
  }
}

float linear_interpolation(float value1, float value2, const float &t,
                           const float &dt, const float &Tprint) {
  return (value2 - value1) * (Tprint - t) / dt + value1;
}

std::string intToStringWithLeadingZeros(int number, int totalLength) {
  return std::string(totalLength - std::to_string(number).length(), '0') +
         std::to_string(number);
}

// //////////////////////////////////////////////////////////////////////////////////////////
void printProgress(double percentage) {
  int val = (int)(percentage * 100);
  int lpad = (int)(percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}