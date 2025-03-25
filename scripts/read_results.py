"""This module contains functions to read the numerical results."""

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import xarray as xr
import os

def read_results(name_run: str) -> xr.Dataset:
    """Read the numerical results.

    Parameters
    ----------
    name_run : str
        Name of the run.

    wave_amplitude : float
        Amplitude of the wave [m].

    wave_period : float
        Period of the wave [s].
    
    Returns
    -------
    ds : xarray.Dataset
        Dataset with the numerical results.
    """
    # folder name for the run
    dir_run = os.path.join('code/results', name_run)

    # read info file binary
    info = np.fromfile(os.path.join(dir_run, 'info.bin'), dtype=np.float32)

    # time arra
    Time = info[0]
    dt = info[1]
    Kprint = int(round(info[2]))
    time = np.linspace(0, Time, Kprint)

    # time array for the gauges
    dt_gauges = info[3]
    K_gauges = int(round(info[4]))
    time_gauges = np.linspace(0, Time, K_gauges)

    path_gauge = os.path.join(dir_run, 'TimeSeries')
    nbr_gauges = len([f for f in os.listdir(path_gauge) if os.path.isfile(os.path.join(path_gauge, f))])

    # x axis
    Lx = info[5]
    dx = info[6]
    Nx = int(round(info[7]))
    x1 = np.linspace(0, Lx, Nx)

    # water depth
    depth = info[8]

    # index of the position of the wave maker
    iwm = int(round(info[9]))

    # change the reference of the x axis
    x1 = x1 - x1[iwm]

    # read the free surface elevation evolution
    eta_all = np.zeros((Kprint, len(x1)))
    for k in range(Kprint):
        eta = np.fromfile(os.path.join(dir_run, 'FreeSurface', f'H_{k:06d}.bin'), dtype=np.float32)
        eta_all[k, :] = eta

    # read the time series at the gauges

    # change the sie of the array if the length of the time series is different -> Rare case due the float point precision
    if nbr_gauges > 0:
        eta = np.fromfile(os.path.join(path_gauge, f'Gauge_{0}.bin'), dtype=np.float32)
        K_gauges = len(eta) - 1
        time_gauges = np.linspace(0, Time, K_gauges)

    eta_gauges = np.zeros((K_gauges, nbr_gauges))
    index_gauges = np.zeros(nbr_gauges)
    for i in range(nbr_gauges):
        eta = np.fromfile(os.path.join(path_gauge, f'Gauge_{i}.bin'), dtype=np.float32)
        index_gauges[i] = round(eta[0])
        eta_gauges[:, i] = eta[1:]

    index_gauges = index_gauges.astype(int)

    # create the xarray dataset
    ds = xr.Dataset(
        {
            'eta': (['time', 'x'], eta_all, {'units': 'm', 'description': 'Free surface elevation'}),
            'eta_gauges': (['time_gauges', 'gauges'], eta_gauges, {'units': 'm', 'description': 'Free surface timeseries at the gauges'}),
            'index_gauges': (['gauges'], index_gauges, {'description': 'Index of the position of the gauges'}),
        },
        coords={
            'time': time,
            'x': x1,
            'time_gauges': time_gauges,
            'gauges': np.arange(nbr_gauges),
        },
        attrs={
            'name_run': name_run,
            'water_depth': depth,
        },
    )

    return ds




    
