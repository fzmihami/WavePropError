#pragma once

/**
 * @file FluxComputation.h
 * @brief Header file defining the flux computation classes and functions.
 *
 * The `flux` class computes the Shallow Water Equations fluxes for a collocated
 * grid using Finite Volume Godunov's method. It includes different
 * reconstruction orders and flux computation schemes. The `flux_staggered`
 * class computes the Shallow Water Equations fluxes for a staggered grid using
 * a conservative staggered scheme (Stelling and Duinmeijer, 2003). This flux
 * computation is based on second-order reconstruction.
 */

#include "TimeVariables.h"

/**
 * @class flux
 * @brief Computes the fluxes for the Shallow Water Equations in a collocated
 * grid.
 *
 * The `flux` class is responsible for computing the fluxes in the shallow water
 * equations using Finite Volume Godunov's method. It includes different
 * reconstruction orders and flux computation schemes. The class also computes
 * the wave speed for the CFL condition.
 */
class flux {
private:
  unsigned int N; ///< Number of grid points

public:
  float *WS;    ///< Wave speed for the CFL computation
  float *Hn;    ///< Water height at the current time step
  float *Qn;    ///< Momentum (Hn*Un) at the current time step
  float *HL;    ///< Water height at the left cell interface
  float *HR;    ///< Water height at the right cell interface
  float *QL;    ///< Momentum (Hn*Un) at the left cell interface
  float *QR;    ///< Momentum (Hn*Un) at the right cell interface
  float *FluxH; ///< Flux of water height at the cell interface
  float *FluxU; ///< Flux of momentum at the cell interface

public:
  /**
   * @brief Constructor for the flux class.
   *
   * Initializes the flux class with the number of grid points and allocates
   * memory for the flux variables.
   *
   * @param N Number of grid points in the computational domain.
   */
  flux(unsigned int N);

  /**
   * @brief Computes the fluxes for the Shallow Water Equations.
   *
   * This function computes the fluxes using the specified reconstruction order
   * and scheme.
   *
   * @param state Reference to a State object containing the current state of
   * the system.
   * @param Scheme The flux computation scheme to be used (e.g., HLL, HLLC).
   * @param OrderReconstruction The order of reconstruction to be used (1, 2, 3,
   * or 5).
   */
  void ComputeFlux(const State &state, const scheme &Scheme,
                   const unsigned int &OrderReconstruction);

  /**
   * @brief Destructor for the flux class.
   *
   * Cleans up dynamically allocated memory for the flux variables.
   */
  ~flux();
};

/**
 * @class flux_staggered
 * @brief Computes the fluxes for the Shallow Water Equations in a staggered
 * grid.
 *
 * The `flux_staggered` class is responsible for computing the fluxes in the
 * shallow water equations using a conservative staggered scheme (Stelling and
 * Duinmeijer, 2003). This flux computation is based on second-order
 * reconstruction.
 */
class flux_staggered {
private:
  unsigned int N; ///< Number of grid points

public:
  float *WS;   ///< Wave speed for the CFL computation
  float *HUn;  ///< Hn*Un approximation in the continuity equation
  float *HUUn; //< Hn*Un approximation in the momentum equation

public:
  /**
   * @brief Constructor for the flux_staggered class.
   *
   * Initializes the flux_staggered class with the number of grid points and
   * allocates memory for the flux variables.
   *
   * @param N Number of grid points in the computational domain.
   */
  flux_staggered(unsigned int N);

  /**
   * @brief Computes the fluxes for the Shallow Water Equations in a staggered
   * grid.
   *
   * This function computes the fluxes using the specified reconstruction order.
   *
   * @param state Reference to a State object containing the current state of
   * the system.
   * @param OrderReconstruction The order of reconstruction to be used (1 or 2).
   *
   * @cite Stelling and Duinmeijer (2003)
   */
  void ComputeFlux(const State &state, const unsigned int &OrderReconstruction);

  /**
   * brief Destructor for the flux_staggered class.
   *
   * Cleans up dynamically allocated memory for the flux variables.
   */
  ~flux_staggered();
};

/**
 * @brief Performs first-order reconstruction of state variables at cell
 * interfaces.
 *
 * This function approximates the state variables at the interface between
 * computational cells using a first-order upwind reconstruction. It assigns
 * values from the left and right neighboring cells, which can be used as inputs
 * for a local Riemann solver in numerical flux computations.
 *
 * The reconstruction follows:
 * \f[
 * U_L^{i+1/2} = U_i
 * \f]
 * \f[
 * U_R^{i+1/2} = U_{i+1}
 * \f]
 * where \( U_L \) and \( U_R \) are the left and right states at the interface.
 *
 * @param[in]  Hn  Pointer to the array of state variables at cell centers (size
 * N).
 * @param[out] Hl  Pointer to the array storing left interface values (size
 * N-1).
 * @param[out] Hr  Pointer to the array storing right interface values (size
 * N-1).
 * @param[in]  N   Number of grid points in the input array Hn.
 *
 * @warning Ensure that `Hl` and `Hr` arrays have been allocated with at least
 * `N-1` elements.
 */
void Recontruction_Order1(float *Hn, float *Hl, float *Hr, unsigned int N);

/**
 * @brief Performs second-order MUSCL reconstruction of state variables at cell
 * interfaces.
 *
 * This function computes the left and right interface states using the
 * Monotonic Upstream-Centered Scheme for Conservation Laws (MUSCL), which
 * provides a second-order accurate reconstruction with a Total Variation
 * Diminishing (TVD) slope limiter.
 *
 * The reconstructed interface values are given by:
 * \f[
 * U_L^{i+1/2} = U_i + \phi (r_L) (U_i - U_{i-1})
 * \f]
 * \f[
 * U_R^{i+1/2} = U_{i+1} - \phi (r_R) (U_{i+2} - U_{i+1})
 * \f]
 *
 * where the slope ratio is defined as:
 * \f[
 * r_L = \frac{U_{i+1} - U_i}{U_i - U_{i-1}}, \quad
 * r_R = \frac{U_{i+1} - U_i}{U_{i+2} - U_{i+1}}
 * \f]
 *
 * The function uses a Total Variation Diminishing (TVD) slope limiter to
 * prevent oscillations. We use the Van-Leer MinMod limiter, defined as: \f[
 * \phi (r, \theta) = \max \left(0, \min \left( \theta r, \frac{1 + r}{2},
 * \theta \right) \right) \f]
 *
 * the value of \theta is deined as a constant in GlobalVariables.cpp file.
 *
 * @param[in]  Hn  Pointer to the array of state variables at cell centers (size
 * N).
 * @param[out] Hl  Pointer to the array storing left interface values (size
 * N-1).
 * @param[out] Hr  Pointer to the array storing right interface values (size
 * N-1).
 * @param[in]  N   Number of grid points in the input array Hn.
 *
 * @warning Ensure `Hl` and `Hr` arrays have at least `N-1` elements to avoid
 * out-of-bounds errors.
 */
void Recontruction_Order2(float *Hn, float *Hl, float *Hr, unsigned int N);

/**
 * @brief Performs third-order TVD reconstruction of state variables at cell
 * interfaces.
 *
 * This function implements the third-order Total Variation Diminishing (TVD)
 * interpolation proposed by Kim and Kim (2005) for numerical flux computation.
 * It ensures stability and prevents spurious oscillations using a nonlinear
 * slope limiter.
 *
 * The reconstructed interface values are given by:
 * \f[
 * U_L^{i+1/2} = U_i + \max \left(0, \min \left(2, 2r_{L,i}, \beta_L \right)
 * \right) (U_i - U_{i-1}) \f] \f[ U_R^{i+1/2} = U_{i+1} - \max \left(0, \min
 * \left(2, 2r_{R,i+1}, \beta_R \right) \right) (U_{i+2} - U_{i+1}) \f]
 *
 * where:
 * \f[
 * \beta_L = \frac{1 + 2r_{L,i}}{3}, \quad
 * \beta_R = \frac{1 + 2r_{R,i+1}}{3}
 * \f]
 *
 * and the slope ratios are:
 * \f[
 * r_{L,i} = \frac{U_{i+1} - U_i}{U_i - U_{i-1}}, \quad
 * r_{R,i+1} = \frac{U_{i+1} - U_i}{U_{i+2} - U_{i+1}}
 * \f]
 *
 * @param[in]  Hn  Pointer to the array of state variables at cell centers (size
 * N).
 * @param[out] Hl  Pointer to the array storing left interface values (size
 * N-1).
 * @param[out] Hr  Pointer to the array storing right interface values (size
 * N-1).
 * @param[in]  N   Number of grid points in the input array Hn.
 *
 * @note Assumes proper boundary handling for the first and last two grid
 * points.
 *
 * @warning Ensure `Hl` and `Hr` arrays have at least `N-1` elements to avoid
 * out-of-bounds errors.
 */
void Recontruction_Order3(float *Hn, float *Hl, float *Hr, unsigned int N);

/**
 * @brief Performs fifth-order TVD reconstruction of state variables at cell
 * interfaces.
 *
 * This function implements a fifth-order Total Variation Diminishing (TVD)
 * interpolation for numerical flux computation. It uses a nonlinear slope
 * limiter to prevent spurious oscillations while maintaining high accuracy.
 *
 * The reconstructed interface values are given by:
 * \f[
 * U_L^{i+1/2} = U_i + \frac{1}{2} \max \left(0, \min \left(2, 2r_{L,i}, \beta_L
 * \right) \right) (U_i - U_{i-1}) \f] \f[ U_R^{i+1/2} = U_{i+1} - \frac{1}{2}
 * \max \left(0, \min \left(2, 2r_{R,i+1}, \beta_R \right) \right) (U_{i+2} -
 * U_{i+1}) \f]
 *
 * where the slope limiter coefficients are:
 * \f[
 * \beta_L = \frac{-2/r_{L,i-1} + 11 + 24r_{L,i} - 3r_{L,i}r_{L,i+1}}{30},
 * \quad \beta_R = \frac{-2/r_{R,i+2} + 11 + 24r_{R,i+1} -
 * 3r_{R,i}r_{R,i+1}}{30} \f]
 *
 * The slope ratios are computed as:
 * \f[
 * r_{L,i} = \frac{U_{i+1} - U_i}{U_i - U_{i-1}},
 * \quad r_{R,i+1} = \frac{U_{i} - U_{i-1}}{U_{i+1} - U_i}
 * \f]
 *
 * @param[in]  Hn  Pointer to the array of state variables at cell centers (size
 * N).
 * @param[out] Hl  Pointer to the array storing left interface values (size
 * N-1).
 * @param[out] Hr  Pointer to the array storing right interface values (size
 * N-1).
 * @param[in]  N   Number of grid points in the input array Hn.
 *
 * @note Assumes proper boundary handling for the first and last two grid
 * points.
 *
 * @warning Ensure `Hl` and `Hr` arrays have at least `N-1` elements to avoid
 * out-of-bounds errors.
 *
 * @cite Kim and Kim (2005)
 */
void Recontruction_Order5(float *Hn, float *Hl, float *Hr, unsigned int N);

/**
 * @brief Computes the MinMod limiter function.
 *
 * The MinMod function is a slope limiter used in high-resolution schemes
 * to prevent spurious oscillations. It selects the smallest absolute
 * value among the inputs if they have the same sign; otherwise, it returns
 * zero.
 *
 * @param a First input value.
 * @param b Second input value.
 * @param c Third input value.
 * @return The MinMod limited value.
 *
 * The function follows these rules:
 * - If all three inputs are positive, it returns the minimum of the three.
 * - If all three inputs are negative, it returns the maximum (most negative) of
 * the three.
 * - Otherwise, it returns zero to avoid introducing new extrema.
 *
 * @cite Kim and Kim (2005)
 */
float minmod(float a, float b, float c);

/**
 * @brief Upwind function for 1D advection problems.
 *
 * This function implements the basic upwind scheme. It returns one of the
 * two input values, depending on the sign of the velocity.
 * - If the velocity (V) is positive, the function returns Q1.
 * - If the velocity (V) is negative, the function returns Q2.
 *
 * @param V The velocity.
 * @param Q1 The value at the left of the interface.
 * @param Q2 The value at the right of the interface.
 * @return The value corresponding to the upwind direction based on the
 * velocity.
 */
float upwind(float V, float Q1, float Q2);

/**
 * @brief Upwind function with slope limiting for 1D advection problems.
 *
 * This function implements an upwind scheme with slope limiting, using a
 * higher-order approximation for the advection of the variable.
 * - If the velocity (V) is positive, the function computes a limited
 *   difference between Q2 and Q3, taking into account the values Q1, Q2,
 *   and Q3.
 * - If the velocity (V) is negative, the function computes a limited
 *   difference between Q3 and Q4, considering Q2, Q3, and Q4.
 *
 * The function uses the MinMod limiter to avoid introducing non-physical
 * oscillations. The limiter is scaled by a factor `thetaRC` which controls the
 * amount of limiting applied. The value of `thetaRC` is defined in the
 * GlobalVariables.cpp file.
 *
 * @param V The velocity.
 * @param Q1 The value at the first neighboring point.
 * @param Q2 The value at the second neighboring point.
 * @param Q3 The value at the third neighboring point.
 * @param Q4 The value at the fourth neighboring point.
 * @return The limited value at the interface.
 */
float upwind(float V, float Q1, float Q2, float Q3, float Q4);

/**
 * @brief Central-upwind scheme for numerical flux computation.
 *
 * This function computes the numerical fluxes for the variables `h` (height)
 and `Q` (momentum)
 * using the central-upwind scheme, as described by Kurganov and Petrova (2007).
 *
 * The fluxes are calculated for each cell interface based on the left and right
 states, `Hl`, `Hr`,
 * `Ql`, and `Qr`, and the speed of sound, `aL` and `aR`. The fluxes are
 computed for both `h` and `Q`
 * using a central-upwind formula with local speeds `aL` and `aR`, and the
 fluxes are stored in
 * `FluxH` and `FluxQ` respectively.
 *
 * @param Hl The array of left state values for the height (h).
 * @param Hr The array of right state values for the height (h).
 * @param Ql The array of left state values for the momentum (Q).
 * @param Qr The array of right state values for the momentum (Q).
 * @param FluxH The array where the numerical flux for height (h) will be
 stored.
 * @param FluxQ The array where the numerical flux for momentum (Q) will be
 stored.
 * @param N The number of cells.
 *

 * The central-upwind numerical flux for each cell interface is computed as:
 *
 * \f[
 * F_{h}^{i+1/2} = \frac{a_R F_{hL} - a_L F_{hR}}{a_R - a_L} + \frac{a_R a_L
 (H_r^{i} - H_l^{i})}{a_R - a_L}
 * \f]
 *
 * \f[
 * F_{Q}^{i+1/2} = \frac{a_R F_{uL} - a_L F_{uR}}{a_R - a_L} + \frac{a_R a_L
 (Q_r^{i} - Q_l^{i})}{a_R - a_L}
 * \f]
 *
 * where:
 *
 * \f[
 * a_R = \max\left(0, u_R + \sqrt{G H_R}, u_L + \sqrt{G H_L}\right)
 * \f]
 *
 * \f[
 * a_L = \min\left(0, u_R - \sqrt{G H_R}, u_L - \sqrt{G H_L}\right)
 * \f]
 *
 * @cite Kurganov, A., & Petrova, G. (2007). Central-upwind schemes for the
 shallow water equations.
 */
void Flux_CentralUpwind(float *Hl, float *Hr, float *Ql, float *Qr,
                        float *FluxH, float *FluxQ, unsigned int N);

/**
 * @brief HLL (Harten, Lax, and van Leer) scheme for numerical flux computation.
 *
 * This function computes the numerical fluxes for the variables `h` (height)
 * and `Q` (momentum) using the HLL scheme (Harten, Lax, and van Leer, 1983).
 *
 * The function calculates the fluxes for each cell interface based on the left
 * and right states and local wave speeds, with the fluxes stored in the arrays
 * `FluxH` and `FluxQ`. The fluxes are computed using the HLL algorithm.
 *
 * The fluxes for height (`h`) and momentum (`Q`) are calculated using the
 * following expressions:
 *
 * \f[
 * F_{h}^{i+1/2} =
 * \begin{cases}
 * F_{hL}, & \text{if } s_L > 0 \\
 * F_{hR}, & \text{if } s_R < 0 \\
 * \frac{s_R F_{hL} - s_L F_{hR} + (s_L s_R) (H_R - H_L)}{s_R - s_L}, &
 * \text{otherwise} \end{cases} \f]
 *
 * \f[
 * F_{Q}^{i+1/2} =
 * \begin{cases}
 * F_{Q_L}, & \text{if } s_L > 0 \\
 * F_{Q_R}, & \text{if } s_R < 0 \\
 * \frac{s_R F_{Q_L} - s_L F_{Q_R} + (s_L s_R) (Q_R - Q_L)}{s_R - s_L}, &
 * \text{otherwise} \end{cases} \f]
 *
 * where:
 * - \( s_L \) is the left wave speed,
 * - \( s_R \) is the right wave speed,
 * - \( F_{hL}, F_{hR} \) are the fluxes for height,
 * - \( F_{Q_L}, F_{Q_R} \) are the fluxes for momentum.
 *
 * The wave speeds \( s_L \) and \( s_R \) are computed based on the local
 * states and the characteristic speed of the system. The value \( h0 \) is
 * computed based on the formulation in Toro (2009).
 *
 * @param Hl The array of left state values for the height (h).
 * @param Hr The array of right state values for the height (h).
 * @param Ql The array of left state values for the momentum (Q).
 * @param Qr The array of right state values for the momentum (Q).
 * @param FluxH The array where the numerical flux for height (h) will be
 * stored.
 * @param FluxQ The array where the numerical flux for momentum (Q) will be
 * stored.
 * @param N The number of cells.
 *
 * @cite Harten, A., Lax, P. D., & van Leer, B. (1983).
 *       "On Upstream Differencing and Godunov-Type Schemes for Hyperbolic
 * Conservation Laws." SIAM Review, 25(1), 35-61.
 */
void Flux_HLL(float *Hl, float *Hr, float *Ql, float *Qr, float *FluxH,
              float *FluxQ, unsigned int N);

/**
 * @brief HLLC (Harten, Lax, van Leer, and Contact) scheme for numerical flux
 * computation.
 *
 * This function computes the numerical fluxes for the variables `h` (height)
 * and `Q` (momentum) using the HLLC scheme (Harten, Lax, van Leer, and Contact,
 * 1983). The HLLC scheme is a modification of the HLL scheme, which
 * incorporates a contact wave and resolves the Riemann problem by considering
 * three waves: a left and right moving wave, and a contact wave.
 *
 * The function computes the fluxes for each cell interface based on the left
 * and right states (`Hl`, `Hr`, `Ql`, `Qr`) and the corresponding wave speeds.
 *
 * The fluxes for height (`h`) and momentum (`Q`) are calculated using the
 * following expressions:
 *
 * \f[
 * F_{h}^{i+1/2} =
 * \begin{cases}
 * F_{hL}, & \text{if } s_L > 0 \\
 * F_{hR}, & \text{if } s_R < 0 \\
 * F_{hL} + s_L (h^*L - h_L), & \text{if } 0 \leq s_0 \\
 * F_{hR} + s_R (h^*R - h_R), & \text{if } s_0 < 0
 * \end{cases}
 * \f]
 *
 * \f[
 * F_{Q}^{i+1/2} =
 * \begin{cases}
 * F_{Q_L}, & \text{if } s_L > 0 \\
 * F_{Q_R}, & \text{if } s_R < 0 \\
 * F_{Q_L} + s_L (q^*L - h_L u_L), & \text{if } 0 \leq s_0 \\
 * F_{Q_R} + s_R (q^*R - h_R u_R), & \text{if } s_0 < 0
 * \end{cases}
 * \f]
 *
 * where:
 * - \( s_L \) is the left wave speed,
 * - \( s_R \) is the right wave speed,
 * - \( s_0 \) is the middle wave speed (contact wave),
 * - \( F_{hL}, F_{hR} \) are the fluxes for height,
 * - \( F_{Q_L}, F_{Q_R} \) are the fluxes for momentum,
 * - \( h^*L, h^*R \) are the star variables for height at the left and right
 * states,
 * - \( q^*L, q^*R \) are the star variables for momentum at the left and right
 * states.
 *
 * The function uses intermediate variables, including \( h0 \) (the averaged
 * height), \( u0 \) (the averaged velocity), and calculates the middle wave
 * speed \( s_0 \) based on the expressions presented in Toro (2001).
 *
 * @param Hl The array of left state values for the height (h).
 * @param Hr The array of right state values for the height (h).
 * @param Ql The array of left state values for the momentum (Q).
 * @param Qr The array of right state values for the momentum (Q).
 * @param FluxH The array where the numerical flux for height (h) will be
 * stored.
 * @param FluxQ The array where the numerical flux for momentum (Q) will be
 * stored.
 * @param N The number of cells.
 *
 * @cite Toro, 2001.
 *       "Shock-capturing methods for free-surface shallow flows."
 */
void Flux_HLLC(float *Hl, float *Hr, float *Ql, float *Qr, float *FluxH,
               float *FluxQ, unsigned int N);
