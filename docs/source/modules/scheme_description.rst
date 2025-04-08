
Description of Numerical Schemes
=================================

In this section, we present the hybrid FV-FD numerical schemes, which will serve as a comparison to the conservative staggered scheme in the solution of wave propagation.\
The numerical solver is designed to solve the 1D Nwogu's equation using various Godunov-type schemes, flux reconstruction methods, and time integration techniques.

The 1D Nwogu equations can be expressed in their conservative form as follows:

.. math::

   \label{Eq_diff:2110}
       U_t + F_x + S = 0

where:

.. math::

   \label{Eq_diff:2111}
       U = 
       \left[\begin{array}{c}
       H \\
       Hu + H\psi_m \end{array} \right] 
       \hspace{20pt}
       F = \left[\begin{array}{c}
       Hu \\
       Hu^2 + \frac{1}{2} gH^2\end{array} \right] 
       \hspace{20pt}
       S = \left[\begin{array}{c}
       \psi_c \\
       -gH h_x +u \psi_c - H_t \psi_m
       \end{array} \right]

.. math::

   \label{Eq_diff:2103}
       \psi_c =
       \Bigg[
       \left( \frac{z_{\alpha}^2}{2} - \frac{h^2}{6} \right) h u_{xx}
       + \left( z_{\alpha} + \frac{h}{2}\right) h \left(hu \right)_{xx}
       \Bigg]_x

.. math::

   \label{Eq_diff:2104}
       \psi_m = \frac{z_{\alpha}^2}{2} u_{xx} +   z_{\alpha} \left( hu \right)_{xx}

The method of characteristics requires the homogeneous part of the equation to be isolated. 
Therefore, to preserve the hyperbolicity of the hydrostatic component, 
the pressure term :math:`gH \eta_x` is split into an flux gradient and a source term.
The Finite Volume Method divides the computational domain into discrete
volumes that store the average values of each conserved variable. These
values are then updated at each time step using the fluxes at the cell
interfaces. The FV method is a cell-centered scheme, meaning that the 
conserved variables are stored at the center of each cell. The conservative
Nwogu's equation is integrated over each cell, and the solution is written as:


.. math::

   \label{Eq_diff:0101}
       \frac{U_i^{n+1} -  U_i^{n}}{\Delta t}
                    +
                   \frac{1}{\Delta x} \left[ F_{i+\frac{1}{2}}^{n} - F_{i-\frac{1}{2}}^{n}
                   \right]
                   + S_{i}^n
                   = 0 .

The source terms in Nwogu's equations :math:`S_{i}^n` are computed with central FD approximations. It's important to note that in the
case of a non-flat bottom, the computation of the bottom variation
term :math:`gHh_x` is not straightforward and requires the use of
specialized techniques to achieve a well-balanced and
positivity-preserving numerical solution (*e.g*., :cite:`audusse2015simple, chertock2015well, liang2009numerical`).




Godunov-type schemes
~~~~~~~~~~~~~~~~~~~~

**Harten, Lax, van Leer (HLL) scheme**

Harten, Lax, and van Leer :cite:`harten1983upstream` introduced a simplified Riemann solver based on the fastest 
and slowest wave speeds, yielding the numerical flux approximation:

.. math::

   \label{Eq_diff:0202}
       F_{i+\frac{1}{2}}^{hll}= 
       \begin{cases}
           F_{L}        & \text{if} \quad 0 \leqslant S_{L}\\
           \frac{S_R F_L - S_L F_R + S_R S_L \left( U_R - U_L \right)}{S_R - S_L}   & \text{if} \quad S_L \leqslant 0 \leqslant S_R\\
           F_{R}        & \text{if} \quad S_R \leqslant 0 
       \end{cases}

The wave speeds :math:`S_L` and :math:`S_R` can be estimated using several approaches; here, we adopted the method proposed by :cite:`toro2001shock`:


.. math::

   \label{Eq_diff:0203}
       S_L = u_L - q_L \sqrt{gH_L} \hspace{10mm}
       S_R = u_R + q_R \sqrt{gH_R}

where:

.. math::

   \label{Eq_diff:0204}
       q_K = 
       \begin{cases}
           \sqrt{\dfrac{1}{2} \left( \dfrac{H_* (H_* + H_K)}{H_K^2} \right)} & \text{if } H_K \leq H_* \\
           1 & \text{if } H_K > H_*
       \end{cases}
       \quad \text{for } K = L, R

and:

.. math::

   \label{Eq_diff:0205}
       H_* = \dfrac{1}{2} (H_L + H_R)
             - \dfrac{1}{4} (H_L + H_R)
             \dfrac{u_R - u_L}{\sqrt{gH_L} + \sqrt{gH_R}}

**Harten, Lax, van Leer Contact (HLLC) scheme**


The HLL scheme was later corrected by :cite:`fraccarollo1995experimental`  to account for the
influence of intermediate waves. The HLLC numerical flux is computed as
follows:

.. math::

   \label{Eq_diff:0102}
       F_{i+\frac{1}{2}}^{hllc}= 
       \begin{cases}
           F_{L}        & \text{if} \quad 0 \leqslant S_{L}\\
           F_{*L} = F_{L} + S_{L} \left(U_{*L} - U_{L}\right)
           & \text{if} \quad S_L \leqslant 0 \leqslant S_{*}\\
           F_{*R} =   F_{R} + S_{R} \left(U_{*R} - U_{R}\right)
           & \text{if} \quad S_{*} \leqslant 0 \leqslant S_R\\
           F_{R}        & \text{if} \quad S_R \leqslant 0 
       \end{cases}

The state :math:`U_{*L}` and :math:`U_{*R}` are given by:

.. math::

   \label{Eq_diff:0103}
       U_{*K} = H_K \left( \frac{S_K - u_K}{S_K - S_{*}} \right)
       \left[\begin{array}{c}
       1 \\
       S^{*}
       \end{array} \right] 
       \hspace{20pt}
       K=L, R

We compute the wave middle wave speed :math:`S_{*}` based on the
formulations provided by :cite:`toro2001shock`:

.. math::

   \label{Eq_diff:0104}
       S_{*} 
       =
       \frac{S_L H_R \left(u_R - S_R \right) - S_R H_L \left(u_L - S_L \right)}{H_R \left(u_R - S_R \right) - H_R \left(u_L - S_L \right)}

**Central-Upwind scheme**

The numerical flux :math:`F_{i+\frac{1}{2}}^{n}` can also be computed
with a Riemann-free solver such as central-upwind scheme (*i.e*.,
:cite:`kurganov2007second`):

.. math::

   \label{Eq_diff:0105}
       F_{i+\frac{1}{2}}^{cu}= 
       \frac{a_{R} F_{L} - a_{L} F_{R}}{a_{R} - a_{L}}
       +
       \frac{a_R a_L}{a_R - a_L}
       \left(
       U_R - U_L
       \right)

where:

.. math::

   \begin{aligned}
   \label{Eq_diff:01055}
       a_{R} &= \max
       \left(
       0, 
       u_R + \sqrt{gH_R},
       u_L + \sqrt{gH_L}
       \right)  \\
       a_{L} &= \min
       \left(
       0, 
       u_R - \sqrt{gH_R},
       u_L - \sqrt{gH_L}
       \right) 
   \end{aligned}

Flux reconstruction
~~~~~~~~~~~~~~~~~~~

To compute the numerical flux, we need to approximate the state
variables at the cell interface. A first-order approximation can be
achieved by taking the values from the left and right cells as input
into the local Riemann solver:

.. math::

   \label{Eq_diff:0106}
       \big[ U_L \big]_{i+\frac{1}{2}} = U_{i} 
       \hspace{30pt}
       \big[ U_R \big]_{i+\frac{1}{2}} = U_{i+1}

For improved accuracy, several high-order schemes have been proposed.
These methods employ flux limiters to achieve oscillation-free
computations:

**Second-order MUSCL reconstruction**

The standard second-order MUSCL reconstruction proposed by :cite:`zhou2001surface` can be described as:

.. math::

   \label{eq_diff:3101}
       \big[ U_L \big]_{i+\frac{1}{2}} = U_i + \frac{1}{2} \phi \left( r_L \right) \left( U_{i} - U_{i-1} \right)

.. math::

   \label{eq_diff:3102}
       \big[ U_R \big]_{i+\frac{1}{2}}  = U_{i+1} - \frac{1}{2} \phi \left( r_R \right) \left( U_{i+2} - U_{i+1} \right)

where:

.. math::

   \label{eq_diff:3103}
       r_L = \frac{U_{i+1} - U_{i}}{U_{i} - U_{i-1}} 
       \hspace{30pt}
       r_R = \frac{U_{i+1} - U_{i}}{U_{i+2} - U_{i+1}}

:math:`\phi` is a TVD slope limiter. Here we use MinMod slope limiter,
which is defined as follows:

.. math::

   \label{eq_diff:3105}
       \phi \left( \boldsymbol{r} \right) = \max \big( 0, \min \left(1 , r \right) \big)


**Third-order MLP reconstruction**

:cite:`kim2005accurate` proposed the multi-dimensional
limiting process (MLP) for multi-dimensional flows. For third-order
accuracy, the flux interpolation is computed as:

.. math::

   \label{eq_diff:3106}
       \big[ U_L \big]_{i+\frac{1}{2}} = U_i + \frac{1}{2} \max \left(0, \min \left(2, 2 r_{L,i}, \beta_L \right) \right) \left( U_{i} - U_{i-1} \right) ,

.. math::

   \label{eq_diff:3107}
       \big[ U_R \big]_{i+\frac{1}{2}} = U_{i+1} - \frac{1}{2} \max \left(0, \min \left(2, 2 r_{R,i+1}, \beta_R \right) \right) \left( U_{i+2} - U_{i+1} \right) .

The values of :math:`\beta_L` and :math:`\beta_R` are calculated as
follows:

.. math::

   \label{eq_diff:3108}
       \beta_L = \frac{1 + 2 r_{L,i}}{3} ,
       \hspace{30pt}
       \beta_R = \frac{1 + 2 r_{R,i+1}}{3} ,

where:

.. math::

   \label{eq_diff:3109}
       r_{L,i} = \frac{U_{i+1} - U_{i}}{U_{i} - U_{i-1}} ,
       \hspace{30pt}
       r_{R,i+1} = \frac{U_{i+1} - U_{i}}{U_{i+2} - U_{i+1}} .

**Fifth-order MLP reconstruction**

A fifth-order reconstruction was also proposed by :cite:`kim2005accurate`:

.. math::

   \label{eq_diff:3110}
       \big[ U_L \big]_{i+\frac{1}{2}} = U_i + \frac{1}{2} \max \left(0, \min \left(2, 2 r_{L,i}, \beta_L \right) \right) \left( U_{i} - U_{i-1} \right)

.. math::

   \label{eq_diff:31101}
       \big[ U_R \big]_{i+\frac{1}{2}} = U_{i+1} - \frac{1}{2} \max \left(0, \min \left(2, 2 r_{R,i+1}, \beta_R \right) \right) \left( U_{i+2} - U_{i+1} \right)

The values of :math:`\beta_L` and :math:`\beta_R` are calculated as
follows:

.. math::

   \label{eq_diff:3111}
       \beta_L = \frac{-2/r_{L, i-1} + 11 + 24 r_{L, i} - 3  r_{L, i}   r_{L, i+1} }{30} 
       \hspace{20pt}
       \beta_R = \frac{-2/r_{R, i+2} + 11 + 24 r_{R, i+1} - 3  r_{R, i}   r_{R, i+1} }{30}

.. math::

   \label{eq_diff:3112}
       r_{L,i} = \frac{U_{i+1} - U_{i}}{U_{i} - U_{i-1}} 
       \hspace{30pt}
       r_{R,i} = \frac{U_{i} - U_{i-1}}{U_{i+1} - U_{i}}

Time-integration
~~~~~~~~~~~~~~~~

The space-discretization described above is coupled with Strong
Stability Preserving (SSP) Runge-Kutta (RK) methods.

The first-order Euler time-integration is computed as:

.. math::

   \label{Eq_diff:0401}
       U_i^{n+1} 
       = U_i^{n}
       - \Delta t
       \left(
       \frac{F_{i+\frac{1}{2}}^{n} - F_{i-\frac{1}{2}}^{n}}{\Delta x}
       + S_i^{n}
       \right)

The second-order RK time integration uses a two-stage time stepping:

.. math::

   \begin{aligned}
   \label{Eq_diff:0402}
       U_i^{(1)} 
       &= U_i^{n}
       - \Delta t
       \left(
       \frac{F_{i+\frac{1}{2}}^{n} - F_{i-\frac{1}{2}}^{n}}{\Delta x}
       + S_i^{n}
       \right) \\
       U_i^{n+1} 
       &= \frac{1}{2}  U_i^{n}
       + \frac{1}{2}  U_i^{(1)}
       - \frac{1}{2} \Delta t
       \left(
       \frac{F_{i+\frac{1}{2}}^{(1)} - F_{i-\frac{1}{2}}^{(1)}}{\Delta x}
       + S_i^{(1)}
       \right)
   \end{aligned}

The third-order RK is expressed as:

.. math::

   \begin{aligned}
   \label{Eq_diff:0403}
       U_i^{(1)} 
       &= U_i^{n}
       - \Delta t
       \left(
       \frac{F_{i+\frac{1}{2}}^{n} - F_{i-\frac{1}{2}}^{n}}{\Delta x}
       + S_i^{n}
       \right) \\
       U_i^{(2)} 
       &= \frac{3}{4}  U_i^{n}
       + \frac{1}{4}  U_i^{(1)}
       - \frac{1}{4} \Delta t
       \left(
       \frac{F_{i+\frac{1}{2}}^{(1)} - F_{i-\frac{1}{2}}^{(1)}}{\Delta x}
       + S_i^{(1)}
       \right) \\
       U_i^{n+1} 
       &= \frac{1}{3}  U_i^{n}
       + \frac{2}{3}  U_i^{(2)}
       - \frac{2}{3} \Delta t
       \left(
       \frac{F_{i+\frac{1}{2}}^{(2)} - F_{i-\frac{1}{2}}^{(2)}}{\Delta x}
       + S_i^{(2)}
       \right) 
   \end{aligned}
