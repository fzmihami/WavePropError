#pragma once
/**
 * @file Wavemaker.h
 * @brief Implements the WaveMaker class for generating wave forcing through a source term.
 *
 * This file provides the implementation of the WaveMaker class, which computes
 * the source term used to generate waves within the computational domain. It also
 * manages and stores the relevant wavemaker parameters. The wave forcing is
 * constructed by superimposing contributions from multiple frequency components.
 * 
 * Reference: Wei et al. (1999)
 */

#include <iostream>

#include "GlobalVariables.h"
#include "ReadUserInput.h"

/**
 * @class WaveMaker
 * @brief Computes the source term for the wavemaker.
 *
 * The `WaveMaker` class is responsible for computing the source term for the
 * wavemaker and adding the term to the momnentum equation.
 */
class WaveMaker {
public:
  unsigned int Nf;    ///< Number of frequencies
  unsigned int N;     ///< Number of grid points
  float Beta;         ///< Beta parameter for the source term
  unsigned int posWM; ///< Position of the wavemaker in the grid
  float depth;        ///< Water depth at the wavemaker
  float width; ///< Width of the wavemaker, computed from the peak wavelength as
               ///< 0.5*Lp
  unsigned int Ns; ///< Number of grid points for the wavemaker

  float Kp; ///< Wavenumber corresponding to the peak frequency

private:
  const float de = 0.5f; ///< Wavemaker constant based on the formulation
                         ///< described by Wei et al. 1999
  const float alpha3 =
      alpha1 +
      float(1.0 /
            3.0); ///<  alpha3 constant defined as alpha1 + 1/3, where alpha1 is
                  ///<  the dispersion coefficient in Nwogu's equation

public:
  float *Freq;  ///< Frequencies as input by the user
  float *Amp;   ///< Amplitudes as input by the user
  float *Phase; ///< Phases as input by the user
  float *Omg;   ///< Angular frequencies computed from the input frequencies
  float *K;     ///< Wavenumbers computed from the angular frequencies
  float *Di;    ///< Magnitude of the source term for each frequency

  float *SourceTerm; ///< Source term for the wavemaker, computed for each grid
                     ///< point inside the wavemaker

public:
  /**
   * @brief Constructor for the WaveMaker class.
   *
   * Initializes the wavemaker parameters and computes the source term for the
   * wavemaker.
   *
   * @param In Reference to an Input object containing user-defined parameters
   * for the wavemaker such as frequencies, amplitudes, and phases.
   * @param N Number of grid points in the computational domain.
   */
  WaveMaker(const Input &In, const unsigned int &N);

  /**
   * @brief Computes the source term that will be added to the momentum
   * equation.
   *
   * This function computes the source term for the wavemaker based on the input
   * parameters and adds it to the momentum equation.
   *
   * @param t Current time step.
   * @param dx Grid spacing in the computational domain.
   */
  void ComputeSourceTerm(const float &t, const float &dx);

  /**
   * @brief Destructor for the WaveMaker class.
   *
   * Cleans up dynamically allocated memory for the wavemaker parameters.
   */
  ~WaveMaker();
};

/**
 * @brief Solves a second-order equation of the form ax^2 + bx + c = 0.
 *
 * This function computes the roots of the quadratic equation using the
 * quadratic formula.
 *
 * @param a Coefficient of x^2
 * @param b Coefficient of x
 * @param c Constant term
 * @return The root of the equation.
 */
float SecondorderEquation(float a, float b, float c);

/**
 * @brief Computes the wavelength from the angular frequency and water depth.
 *
 * This function calculates the wavelength based on the dispersion relation from
 * Nwogu's equation.
 *
 * @param omega Angular frequency
 * @param depth Water depth
 * @return The computed wavelength.
 */
float wavelength(float omega, float depth);