"""
This module contains functions to create empirical spectral distributions from user input data.

It defines an enumeration of common wave spectral models and provides functions to generate
spectral densities, amplitude distributions, and phase distributions for ocean waves.

Functions
---------
empirical_spec :
    Returns the spectral density of a given wave spectrum input.

create_empirical_spec :
    Creates the empirical spectrum and optionally plots it.

Classes
-------
spec : Enum
    Defines the types of wave spectra: Pierson-Moskowitz, JONSWAP, and TMA.

Constants
---------
G : float
    Gravitational acceleration [m/s^2].
"""

import numpy as np
from enum import Enum
import os
import matplotlib.pyplot as plt


G = 9.81  # m/s^2


# Define an enumeration class for wave spectral models
class spec(Enum):
    pierson_moskowitz = 1  # Pierson-Moskowitz spectrum (index 1)
    jonswap = 2  # JONSWAP spectrum (index 2)
    tma = 3  # TMA spectrum (index 3)


# Retrieve specific spectral models
pierson_moskowitz = (
    spec.pierson_moskowitz
)  # Assign Pierson-Moskowitz spectrum to the variable
jonswap = spec.jonswap  # Assign JONSWAP spectrum to the variable
tma = spec.tma  # Assign TMA spectrum to the variable


def empirical_spec(f, fp, spec_type, Hs, water_depth):
    """Returns the spectral density of a given wave spectrum input.

    Parameters
    ----------
    f : array_like
        Frequency vector [Hz].
    fp : float
        Peak frequency [Hz].
    spec_type : spec
        Type of wave spectrum.
    Hs : float
        Significant wave height [m].
    water_depth : float
        Water depth [m].

    Returns
    -------
    S : array_like
        Spectral density [m^2/(Hz)].
    """
    # define the constants
    alpha = 0.0081
    gamma = 3.3
    sigma = np.where(f <= fp, 0.07, 0.09)
    d = water_depth
    Xa = (2 * np.pi * f) * np.sqrt(d / G)
    phi = np.where(Xa < 1, 0.5 * Xa**2, np.where(Xa >= 2, 1, 1 - 0.5 * (2 - Xa) ** 2))

    # define the spectra density
    S0 = (
        alpha
        * G**2
        * (2 * np.pi) ** (-4)
        * f ** (-5)
        * np.exp(-5.0 / 4 * (fp / f) ** 4)
    )

    if spec_type == pierson_moskowitz:
        S = S0

    elif spec_type == jonswap:
        S = S0 * gamma ** (np.exp(-0.5 * ((f - fp) / (sigma * fp)) ** 2))

    elif spec_type == tma:
        S = S0 * gamma ** (np.exp(-0.5 * ((f - fp) / (sigma * fp)) ** 2)) * phi

    else:
        raise ValueError("Invalid spectrum type.")

    # Correct spectrum to get the desired Hs
    S = S * Hs**2 / (16 * np.trapz(S, f))

    return S


def create_empirical_spec(
    spec_type, Hs, Tp, water_depth, run_time, show_plot=False
) -> np.ndarray:
    """create the text file containing the empirical spectrum, which is used as input for the wavemaker.

    Parameters
    ----------
    spec_type : spec
        Type of empirical spectrum.
    Hs : float
        Significant wave height [m].
    Tp : float
        Peak period [s].
    water_depth : float
        Water depth [m].
    run_time : float
        Simulation run time [s].
    show_plot : bool, optional
        Show the plot of the empirical spectrum.

    Returns
    -------
    spec_array : array_like [Nf, 4]
        Empirical spectrum [frequency, density, amplitude, phase].
        units: [Hz, m^2/Hz, m, rad].
    """
    # define the frequency vector
    Lmin = 2.0 * water_depth
    omg_max = np.sqrt(
        G * 2.0 * np.pi / Lmin * np.tanh(2.0 * np.pi / Lmin * water_depth)
    )
    fmax = omg_max / (2.0 * np.pi)

    fmin = max(1.0 / (2.0 * Tp), 1.0 / 30.0)

    Nf = int(np.ceil(run_time * (fmax - fmin))) + 1

    f = np.linspace(fmin, fmax, Nf)

    df = f[1] - f[0]

    # define the peak frequency
    fp = 1.0 / Tp

    # create the empirical spectrum
    S = empirical_spec(f, fp, spec_type, Hs, water_depth)

    # compute the wave amplitude
    amp = np.sqrt(2 * S * df)

    # compute the wave phase
    np.random.seed(0)  # repeatable results
    phase = np.random.uniform(0, 2 * np.pi, len(f))

    # plot spectral density, amplitude and phase
    fig, axs = plt.subplots(3, 1, figsize=(8, 10))
    axs[0].plot(f, S)
    axs[0].set_ylabel("Spectral density [m^2/Hz]")
    axs[0].set_title("Spectral density distribution")
    axs[0].grid(True)
    plt.setp(axs[0].get_xticklabels(), visible=False)

    axs[1].plot(f, amp)
    axs[1].set_ylabel("Amplitude [m]")
    axs[1].set_title("Amplitude distribution")
    axs[1].grid(True)
    plt.setp(axs[0].get_xticklabels(), visible=False)

    axs[2].plot(f, phase)
    axs[2].set_xlabel("Frequency [Hz]")
    axs[2].set_ylabel("Phase [rad]")
    axs[2].set_title("Phase distribution")
    axs[2].grid(True)

    plt.tight_layout()

    if not show_plot:
        plt.close()

    # create the empirical spectrum array
    spec_array = np.column_stack((f, S, amp, phase))

    return spec_array
