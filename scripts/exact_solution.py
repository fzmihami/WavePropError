"""
This module contains functions to compute the exact solution for wave propagation in a channel based on linearized Nwogu's equations, as well as the Airy wave theory.

It includes functions to compute both the exact solution for a single sine wave and spectral wave propagation.

Constants
---------
G : float
    Acceleration due to gravity [m/s^2].
alpha1 : float
    Constant used in Nwogu's equations to achieve optimized dispersion properties for kh in the range [0, pi].
    default: alpha1 = -0.39

Functions
---------
sine_wave_exact_solution :
    Computes the exact solution for a sine wave based on linearized Nwogu's equations or Airy wave theory.

spectral_wave_exact_solution :
    Computes the exact solution for spectral wave propagation using either Nwogu's equations or Airy wave theory.

wavelength_nwogu :
    Computes the wavelength based on the dispersion relation of linearized Nwogu's equations.

waveperiod_nwogu :
    Computes the wave period based on the dispersion relation of linearized Nwogu's equations.

wavelength_airy :
    Computes the wavelength based on the exact Airy wave theory.

Classes
-------
exact_type : Enum
    Enum to choose between Airy and Nwogu's exact solution types.
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scripts.empirical_spec import *

G = 9.81  # m/s^2
alpha1 = -0.39  # Nwogu's constant


class exact_type(Enum):
    airy = 1
    nwogu = 2


def sine_wave_exact_solution(
    water_depth: float,
    wave_amplitude: float,
    wave_period: float,
    time: np.ndarray,
    x: np.ndarray,
    time_gauges: np.ndarray,
    index_gauges: np.ndarray,
    exact_type: exact_type = exact_type.nwogu,
) -> xr.Dataset:
    """Compute the exact solution for wave propagation in a channel based on linearized Nwogu's equations.

    Parameters
    ----------
    water_depth : float
        Water depth [m].
    wave_amplitude : float
        Wave amplitude [m].
    wave_period : float
        Wave period [s].
    time : np.ndarray
        Time array [s].
    x : np.ndarray
        x axis [m].
    time_gauges : np.ndarray
        Time array for the gauges [s].
    index_gauges : np.ndarray
        Locations of the gauges as indices.
    exact_type : exact_type
        Type of exact solution to use.
        default: exact_type.nwogu

    Returns
    -------
    ds : xarray.Dataset
        Dataset with the exact solution.
    """
    # compute the wavelength
    if exact_type == exact_type.airy:
        wavelength = wavelength_airy(wave_period, water_depth)
    elif exact_type == exact_type.nwogu:
        wavelength = wavelength_nwogu(wave_period, water_depth)
    else:
        raise ValueError("Invalid exact solution type")

    # compute the wave number
    k = 2 * np.pi / wavelength

    # compute the angular frequency
    omega = 2 * np.pi / wave_period

    # free surface elevation
    eta_all = np.zeros((len(time), len(x)))
    for i in range(len(time)):
        eta_all[i, :] = wave_amplitude * np.sin(k * x - omega * time[i])

    # timeseries at the gauges
    eta_gauges = np.zeros((len(time_gauges), len(index_gauges)))
    for i in range(len(time_gauges)):
        eta_gauges[i, :] = wave_amplitude * np.sin(
            k * x[index_gauges] - omega * time_gauges[i]
        )

    # create the dataset
    ds = xr.Dataset(
        {
            "eta": (
                ["time", "x"],
                eta_all,
                {"units": "m", "description": "Free surface elevation"},
            ),
            "eta_gauges": (
                ["time_gauges", "gauges"],
                eta_gauges,
                {"units": "m", "description": "Free surface timeseries at the gauges"},
            ),
            "index_gauges": (
                ["gauges"],
                index_gauges,
                {"description": "Index of the position of the gauges"},
            ),
        },
        coords={
            "time": time,
            "x": x,
            "time_gauges": time_gauges,
            "gauges": np.arange(len(index_gauges)),
        },
        attrs={
            "name_run": "exact_solution",
            "water_depth": water_depth,
            "wave_amplitude": wave_amplitude,
            "wave_period": wave_period,
        },
    )

    return ds


def spectral_wave_exact_solution(
    water_depth: float,
    spec_type: spec,
    significant_wave_height: float,
    peak_wave_period: float,
    time: np.ndarray,
    x: np.ndarray,
    time_gauges: np.ndarray,
    index_gauges: np.ndarray,
    exact_type: exact_type = exact_type.nwogu,
) -> xr.Dataset:
    """Compute the exact solution for spectral wave propagation in a channel based on linearized Nwogu's equations.

    Parameters
    ----------
    water_depth : float
        Water depth [m].
    spec_type : spec
        Type of empirical spectrum.
    significant_wave_height : float
        Significant wave height [m].
    peak_wave_period : float
        Peak wave period [s].
    time : np.ndarray
        Time array [s].
    x : np.ndarray
        x axis [m].
    time_gauges : np.ndarray
        Time array for the gauges [s].
    index_gauges : np.ndarray
        Locations of the gauges as indices.
    exact_type : exact_type
        Type of exact solution to use.
        default: exact_type.nwogu

    Returns
    -------
    ds : xarray.Dataset
        Dataset with the exact solution.
    """

    # compute the empirical spectrum
    spec_array = create_empirical_spec(
        spec_type, significant_wave_height, peak_wave_period, water_depth, time[-1]
    )
    freq = spec_array[:, 0]
    density = spec_array[:, 1]
    amp = spec_array[:, 2]
    phase = spec_array[:, 3]

    k = np.zeros_like(freq)
    for i in range(len(freq)):

        if exact_type == exact_type.airy:
            wavelength = wavelength_airy(1.0 / freq[i], water_depth)
        elif exact_type == exact_type.nwogu:
            wavelength = wavelength_nwogu(1.0 / freq[i], water_depth)
        else:
            raise ValueError("Invalid exact solution type")

        k[i] = 2 * np.pi / wavelength

    # compute the free surface elevation
    eta_all = np.zeros((len(time), len(x)))
    for i in range(len(time)):
        for j in range(len(x)):
            eta_all[i, j] = np.sum(
                amp * np.sin(k * x[j] - 2 * np.pi * freq * time[i] + phase)
            )

    # compute the timeseries at the gauges
    eta_gauges = np.zeros((len(time_gauges), len(index_gauges)))
    for i in range(len(time_gauges)):
        for j in range(len(index_gauges)):
            eta_gauges[i, j] = np.sum(
                amp
                * np.sin(
                    k * x[index_gauges[j]] - 2 * np.pi * freq * time_gauges[i] + phase
                )
            )

    # create the dataset
    ds = xr.Dataset(
        {
            "eta": (
                ["time", "x"],
                eta_all,
                {"units": "m", "description": "Free surface elevation"},
            ),
            "eta_gauges": (
                ["time_gauges", "gauges"],
                eta_gauges,
                {"units": "m", "description": "Free surface timeseries at the gauges"},
            ),
            "index_gauges": (
                ["gauges"],
                index_gauges,
                {"description": "Index of the position of the gauges"},
            ),
            "spectral_density": (
                ["frequency"],
                density,
                {"units": "m^2/Hz", "description": "Spectral density"},
            ),
        },
        coords={
            "time": time,
            "x": x,
            "time_gauges": time_gauges,
            "gauges": np.arange(len(index_gauges)),
            "frequency": freq,
        },
        attrs={
            "name_run": "exact_solution",
            "water_depth": water_depth,
            "significant_wave_height": significant_wave_height,
            "peak_wave_period": peak_wave_period,
        },
    )

    return ds


def wavelength_nwogu(wave_period: float, water_depth: float) -> float:
    """Compute the wavelength based on the dispersion relation of linearized Nwogu's equations.

    Parameters
    ----------
    wave_period : float
        Wave period [s].
    water_depth : float
        Water depth [m].

    Returns
    -------
    wavelength : float
        Wavenumber [m].
    """
    omega = 2 * np.pi / wave_period
    beta1 = alpha1 + 1.0 / 3.0

    # solve second order equation with a b and c coefficients
    a = G * beta1
    b = -(water_depth * omega**2 * alpha1 + G)
    c = water_depth * omega**2
    delta = b**2 - 4 * a * c

    # compute the wavenumber
    k = (-b - np.sqrt(delta)) / (2.0 * a)
    k = np.sqrt(k) / water_depth

    # compute the wavelength
    return 2 * np.pi / k


def waveperiod_nwogu(wavelength: float, water_depth: float) -> float:
    """Compute the wave period based on the dispersion relation of linearized Nwogu's equations.

    Parameters
    ----------
    wavelength : float
        Wavelength [m].
    water_depth : float
        Water depth [m].

    Returns
    -------
    wave_period : float
        Wave period [s].
    """
    d = water_depth
    beta1 = alpha1 + 1.0 / 3
    K1 = 2 * np.pi / wavelength
    KH2 = (K1 * d) ** 2

    omega = G * K1**2 * d * (1 - beta1 * KH2) / (1 - alpha1 * KH2)

    omega = np.sqrt(omega)

    return 2 * np.pi / omega


def wavelength_airy(
    wave_period: float, water_depth: float, tol=1e-6, max_iter=1000
) -> float:
    """Compute the wavelength based on the exact Airy wave theory.

    Parameters
    ----------
    wave_period : float
        Wave period [s].
    water_depth : float
        Water depth [m].
    tol : float, optional
        Tolerance for the iteration.
        default: 1e-6
    max_iter : int, optional
        Maximum number of iterations.
        default: 1000

    Returns
    -------
    wavelength : float
        Wavenumber [m].
    """
    omega = 2 * np.pi / wave_period
    k = omega**2 / G  # Initial guess for k

    # Iterate to refine the value of k using the dispersion relation
    for _ in range(max_iter):
        k_new = omega**2 / (G * np.tanh(k * water_depth))
        if np.abs(k - k_new) < tol:
            k = k_new
            break
        k = k_new

    lambda_ = 2 * np.pi / k
    return lambda_
