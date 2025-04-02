"""
This module provides functions for computing errors in the amplitude and phase of the numerical scheme.
The errors are calculated by comparing numerical wave solutions against exact solutions.
The module includes a method for calculating the time shift between two time series based on cross-correlation
and computing the amplitude of a time series using the Hilbert transform.

Functions
---------
compute_time_shift :
    Computes the time shift between two time series using cross-correlation.

compute_amplitude :
    Computes the mean amplitude of a time series using the Hilbert transform.

diffusive_rmse :
    Computes the diffusive RMSE error for the numerical solution.

dispersive_rmse :
    Computes the dispersive RMSE error for the numerical solution.

"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scripts.empirical_spec import *

from scipy.signal import hilbert, detrend
from scipy.signal import correlate


def compute_time_shift(ts1, ts2, time):
    """Compute the time shift between two time series.

    Parameters
    ----------
    ts1 : array_like
        Time series 1.

    ts2 : array_like
        Time series 2.

    time : array_like
        Time array.

    Returns
    -------
    time_shift : float
        Time shift between the two time series.
    """

    # Compute the cross-correlation between the two time series
    xcorr = correlate(ts1, ts2, mode="full")

    # Find the index of the maximum correlation
    max_index = np.argmax(xcorr)

    # Compute the time shift
    time_shift = max_index - (len(ts1) - 1)

    dt = time[1] - time[0]

    return time_shift * dt


def compute_amplitude(ts1):
    """Compute the amplitude of a time series.

    Parameters
    ----------
    ts1 : array_like
        Time series.

    Returns
    -------
    amplitude : float
        Amplitude of the time series.
    """

    # Compute the Hilbert transform of the time series
    analytic_signal = hilbert(ts1)

    # Compute the amplitude of the time series
    amplitude = np.abs(analytic_signal)

    # Compute the mean amplitude
    amplitude = np.mean(amplitude)

    return amplitude


def diffusive_rmse(ds: xr.Dataset):
    """Compute the diffusive rmse error of the numerical solution.

    Parameters
    ----------
    ds_list : array_like [xarray.Dataset]
        Dataset with the numerical results.

    Returns
    -------
    error : float
        Diffusive RMSE error of the numerical solution, calculated by considering the error for all gauges and scaling by the initial amplitude.
    """

    amp = np.zeros(len(ds.gauges))
    amp0 = ds.wave_amplitude

    for i, gauge in enumerate(ds.gauges):
        eta = ds.eta_gauges[:, i]

        # Compute the amplitude of the numerical solution
        amp[i] = compute_amplitude(eta)

    # Compute the diffusive error RMSE
    rmse = np.sqrt(np.mean((amp - amp0) ** 2))
    rmse = rmse / amp0  # scaling by the initial amplitude

    return rmse


def dispersive_rmse(ds: xr.Dataset, ds_exact: xr.Dataset):
    """Compute the dispersive rmse error of the numerical solution.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with the numerical results.

    ds_exact : xarray.Dataset
        Dataset with the exact results.

    Returns
    -------
    error : float
        Dispersive RMSE error of the numerical solution, calculated by considering the error for all gauges and scaling by the wave period.
    """

    time_shift = np.zeros(len(ds.gauges))
    period = ds.wave_period

    for i, gauge in enumerate(ds.gauges):
        eta = ds.eta_gauges[:, i]
        eta_exact = ds_exact.eta_gauges[:, i]

        # Compute the time shift between the numerical and exact solutions
        time_shift[i] = compute_time_shift(eta_exact, eta, ds.time_gauges)

        # correct the time shift
        if time_shift[i] > 0.5 * period:
            time_shift[i] = time_shift[i] - period
        elif time_shift[i] < -0.5 * period:
            time_shift[i] = time_shift[i] + period

        if time_shift[i] > 0.5 * period or time_shift[i] < -0.5 * period:
            time_shift[i] = 0.5 * period * np.sign(time_shift[i])

    # Compute the dispersive error RMSE
    rmse = np.sqrt(np.mean((time_shift - 0) ** 2))
    rmse = rmse / period  # scaling by the wave period

    return rmse
