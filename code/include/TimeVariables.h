#pragma once

/**
 * @file TimeVariables.h
 * @brief Header file defining the State class and functions for managing time
 * variables in the numerical wave solver.
 *
 * The State class represents the state variables (i.e., water height, velocity,
 * etc.) at the current time step. It includes methods for updating the state
 * and managing memory. The file also contains functions for printing progress
 * the progress of the simulation, recording time series data at defined gauges,
 * and saving the simulation results to binary files.
 */

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include "GlobalVariables.h"
#include "ReadUserInput.h"

/**
 * @class State
 * @brief Represents the state variables of the numerical wave solver.
 *
 * The `State` class contains the state variables (i.e., water height, velocity,
 * etc.) at the current time step. It includes methods for updating the state
 * and managing memory.
 */
class State {
public:
  unsigned int N; ///< Number of grid points
  float dx;       ///< Grid spacing
  float dx_1;     ///< Inverse of grid spacing
  float depth;    ///< Water depth, which is constant in the domain
  grid Grid;      ///< Grid type (collocated or staggered)

public:
  float *Hn;  ///< Water height at the current time step
  float *Eta; ///< Water surface elevation at the current time step
  float *Zn;  ///< Bottom elevation at the current time step
  float
      *Pn; ///< Conserved momentum in Nwogu's equation at the current time step
  float *Un; ///< Velocity at the current time step

public:
  /**
   * @brief Constructor for the State class.
   *
   * Initializes the state variables using the initial conditions provided in
   * the input file.
   *
   * @param In Reference to an Input object containing user-defined parameters
   * for the simulation.
   */
  State(const Input &In);

  /**
   * @brief Destructor for the State class.
   *
   * Cleans up dynamically allocated memory for the state variables.
   */
  ~State();

  /**
   * @brief Swaps the pointers of two State objects.
   *
   * This function swaps the pointers of the state variables between two State
   * objects, effectively updating their values. This is useful for updating the
   * state variables at each time step without copying the data.
   *
   * @param other Reference to another State object to swap with.
   */
  void Update(State &other);
};

/**
 * @brief Saves the simulation results (Water surface elevation) to binary files
 * at the current time step.
 *
 * This function saves the water surface elevation data to binary files at the
 * current time step. The files are named according to the current time step and
 * are stored in the specified output folder.
 *
 * @param S0 Reference to the State object containing the current state
 * variables.
 * @param S2 Reference to the State object containing the newly computed state
 * variables: \f$ t + \Delta t \f$.
 * @param In Reference to the Input object containing user-defined parameters
 * for the simulation.
 * @param t1 Current time step.
 * @param dt1 Time step size.
 * @param Kprint Counter for the number of printouts.
 * @param Tprint Time interval for printouts.
 *
 * @note The function performs linear interpolation between the current and new
 * state variables to obtain the water surface elevation at the specified output
 * time.
 */
void PrintVariables(State &S0, State &S2, const Input &In, const float &t1,
                    const float &dt1, unsigned int &Kprint, float &Tprint);

/**
 * @brief Prints the simulation parameters to a binary file.
 *
 * This function prints the simulation parameters to a binary file for reference
 * and documentation purposes. The parameters include the grid size, time step,
 * and other relevant information.
 *
 * @param In Reference to the Input object containing user-defined parameters
 * for the simulation.
 */
void printINFO(const Input &In);

/**
 * @brief Records time series data at specified gauges.
 *
 * This function records the water surface elevation at specified gauges over
 * time. The data is stored in a vector for each gauge and can be printed later.
 *
 * @param S0 Reference to the State object containing the current state
 * variables.
 * @param S2 Reference to the State object containing the newly computed state
 * variables: \f$ t + \Delta t \f$.
 * @param In Reference to the Input object containing user-defined parameters
 * for the simulation.
 * @param Hgauges Pointer to a vector of floats for storing the water surface
 * elevation at the gauges.
 * @param t1 Current time step.
 * @param dt1 Time step size.
 * @param Kgauges Counter for the number of gauges.
 * @param Tgauges Time interval for recording gauges.
 *
 * @note The function performs linear interpolation between the current and new
 * state variables to obtain the water surface elevation at the specified
 * gauges.
 */
void RecordTimeSeries(State &S0, State &S2, const Input &In,
                      std::vector<float> *Hgauges, const float &t1,
                      const float &dt1, unsigned int &Kgauges, float &Tgauges);

/**
 * @brief Prints the recorded time series data to binary files.
 *
 * This function prints the recorded time series data for each gauge to binary
 * files. The files are named according to the gauge index and are stored in the
 * specified output folder.
 *
 * @param Hgauges Pointer to a vector of floats for storing the water surface
 * elevation at the gauges.
 * @param In Reference to the Input object containing user-defined parameters
 * for the simulation.
 */
void PrintTimeSeries(std::vector<float> *Hgauges, const Input &In);

/**
 * @brief Performs linear interpolation between two state variables.
 *
 * This function computes the interpolated value of a state variable at a given
 * time step using linear interpolation.
 *
 * @param H1 Pointer to the state variable at the current time step. It will be
 * updated with the interpolated value.
 * @param H2 Pointer to the state variable at the next time step: \f$ t + \Delta
 * t \f$.
 * @param N Number of grid points.
 * @param t Current time step.
 * @param dt Time step size.
 * @param Tprint Time interval for printouts specified in the input file.
 */
void linear_interpolation(float *H1, float *H2, unsigned int N, const float &t,
                          const float &dt, const float &Tprint);

/**
 * @brief Performs linear interpolation between two float values.
 *
 * This function computes the interpolated value of a float variable at a given
 * time step using linear interpolation.
 *
 * @param value1 The value at the current time step.
 * @param value2 The value at the next time step: \f$ t + \Delta t \f$.
 * @param t Current time step.
 * @param dt Time step size.
 * @param Tprint Time interval for printouts specified in the input file.
 * @return The interpolated value at the specified time.
 */
float linear_interpolation(float value1, float value2, const float &t,
                           const float &dt, const float &Tprint);

/**
 * @brief Converts an integer to a string with leading zeros.
 *
 * This function converts an integer to a string and pads it with leading zeros
 * to ensure it has a specified total length.
 *
 * @param number The integer to convert.
 * @param totalLength The total length of the resulting string (including
 * leading zeros).
 * @return The formatted string with leading zeros.
 */
std::string intToStringWithLeadingZeros(int number, int totalLength);

/**
 * @brief Prints the progress of the simulation to the console.
 *
 * This function displays a progress bar in the console, indicating the
 * percentage of completion of the simulation.
 *
 * @param percentage The percentage of completion (0.0 to 1.0).
 */
void printProgress(double percentage);