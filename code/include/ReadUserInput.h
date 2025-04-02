#pragma once

/**
 * @file ReadUserInput.h
 * @brief Header file defining the Input class and functions for reading user
 * input in the numerical wave solver.
 *
 * The Input class contains methods for reading user-defined parameters from an
 * input file, including simulation parameters, grid properties, wave data, and
 * scheme properties.
 */

#include <cstring>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "GlobalVariables.h"

/**
 * @class Input
 * @brief Represents the user-defined input parameters for the numerical wave
 * solver.
 *
 * The `Input` class contains methods for reading user-defined parameters from
 * an input file, including simulation parameters, grid properties, wave data,
 * and scheme properties. It also includes a method for logging the input
 * parameters to the console.
 */
class Input {
public:
  scheme Scheme; ///< Numerical scheme used for the simulation
  unsigned int
      OrderReconstruction; ///< Order of reconstruction for the numerical scheme
  unsigned int
      OrderTime;       ///< Order of time integration for the numerical scheme
  float CourantNumber; ///< Courant number which will be used to compute the
                       ///< time step
  float dx;            ///< Grid spacing
  float LengthDomain;  ///< Length of the domain
  float WaterDepth;    ///< Water depth, which is constant in the domain
  float Time;          ///< Total computation time
  float dtOutput;      ///< Time step for the water surface elevation output
  float dtGauges;      ///< Time step for the gauges
  std::vector<unsigned int> IndexGauges; ///< Vector of gauges index
  std::vector<std::vector<float>>
      WaveInput;         ///< Vector of wave data (frequency, amplitude, phase)
  std::string ResFolder; ///< Output folder for the simulation results

public:
  /**
   * @brief Constructor for the Input class.
   *
   * Initializes the input parameters by reading them from a file.
   *
   * @param nametest Name of the input file containing user-defined parameters.
   */
  Input(std::string nametest);

  /**
   * @brief Logs the input parameters to the console.
   *
   * This method prints the user-defined parameters to the console for reference
   * and debugging purposes.
   */
  void LOGOUT();
};

/**
 * @brief Converts a string to an enumeration value of the scheme type.
 *
 * This function takes a string representation of a numerical scheme and
 * converts it to the corresponding enumeration value.
 *
 * @param string1 The string representation of the numerical scheme.
 * @return The corresponding enumeration value of the scheme type.
 */
scheme convertScheme(std::string string1);

/**
 * @brief Converts a string to a vector of unsigned integers representing gauge
 * indices.
 *
 * This function takes a string representation of gauge indices and converts it
 * to a vector of unsigned integers.
 *
 * @param string1 The string representation of the gauge indices.
 * @param dx The grid spacing used for the simulation.
 * @return A vector of unsigned integers representing the gauge indices.
 */
std::vector<unsigned int> convertIndexGauges(std::string string1, float dx);
