
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
    xcorr = correlate(ts1, ts2, mode='full')

    # Find the index of the maximum correlation
    max_index = np.argmax(xcorr)
    
    # Compute the time shift
    time_shift = max_index - (len(ts1) - 1)

    dt = time[1] - time[0]

    return time_shift*dt



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


def diffusive_rmse(ds:xr.Dataset):
    """Compute the diffusive rmse error of the numerical solution.
    
    Parameters
    ----------
    ds_list : array_like [xarray.Dataset]
        Dataset with the numerical results.
    
    Returns
    -------
    error : array_like
        Diffusive RMSE error of the numerical solution.
    """
    
    amp = np.zeros(len(ds.gauges))
    amp0 = ds.wave_amplitude

    for i, gauge in enumerate(ds.gauges):
        eta = ds.eta_gauges[:, i]
        
        # Compute the amplitude of the numerical solution
        amp[i] = compute_amplitude(eta)
        
    # Compute the diffusive error RMSE
    rmse = np.sqrt(np.mean((amp - amp0)**2))
    rmse = rmse / amp0 #scaling by the initial amplitude

    return rmse


def dispersive_rmse(ds:xr.Dataset, ds_exact:xr.Dataset):
    """Compute the dispersive rmse error of the numerical solution.
    
    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with the numerical results.
    
    ds_exact : xarray.Dataset
        Dataset with the exact results.
    
    Returns
    -------
    error : array_like
        Dispersive RMSE error of the numerical solution.
    """
    
    time_shift = np.zeros(len(ds.gauges))
    period = ds.wave_period
    
    for i, gauge in enumerate(ds.gauges):
        eta = ds.eta_gauges[:, i]
        eta_exact = ds_exact.eta_gauges[:, i]
        
        # Compute the time shift between the numerical and exact solutions
        time_shift[i] = compute_time_shift(eta_exact, eta, ds.time_gauges)

        # correct the time shift
        if time_shift[i] > 0.5*period:
            time_shift[i] = time_shift[i] - period
        elif time_shift[i] < -0.5*period:
            time_shift[i] = time_shift[i] + period

        if time_shift[i] > 0.5*period or time_shift[i] < -0.5*period:
            time_shift[i] = 0.5*period*np.sign(time_shift[i])
        
    # Compute the dispersive error RMSE
    rmse = np.sqrt(np.mean((time_shift - 0)**2))
    rmse = rmse / period #scaling by the wave period

    return rmse