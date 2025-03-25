"""This module contains functions to run the wave numerical model."""

import numpy as np
import matplotlib.pyplot as plt
import os
from enum import Enum
from scripts.empirical_spec import *
from scripts.read_results import *
import subprocess

G = 9.81  # m/s^2


class scheme(Enum):
    conservative_staggered = 1
    hllc = 2
    central_upwind = 3
    hll = 4


conservative_staggered = scheme.conservative_staggered
hllc = scheme.hllc
central_upwind = scheme.central_upwind
hll = scheme.hll


def run_sine_wave(
    name_run: str,
    scheme: scheme,
    order_reconstruction: int,
    order_time_integration: int,
    courant_number: float,
    grid_size: float,
    domain_size: float,
    water_depth: float,
    wave_amplitude: float,
    wave_period: float,
    run_time: float,
    output_interval: float,
    gauges_locations: np.ndarray,
    gauges_dt: float,
) -> xr.Dataset:
    """Run the wave model with a sine wave input and return the numerical results as a xarray.Dataset.

    Parameters
    ----------
    name_run : str
        Name of the run.
    scheme : scheme
        Numerical scheme used in the model.
    order_reconstruction : int
        Order of the flux reconstruction.
    order_time_integration : int
        Order of the time integration.
    courant_number : float
        Courant number.
    grid_size : float
        Grid size [m].
    domain_size : float
        Domain size [m].
    water_depth : float
        Water depth [m].
    wave_amplitude : float
        Amplitude of the wave [m].
    wave_period : float
        Period of the wave [s].
    run_time : float
        Simulation run time [s].
    output_interval : float
        Output interval [s].
    gauges_locations : array_like [Ng]
        Locations of the gauges [m]. 
        The reference point is the center of the wavemaker.
        Ng is the number of gauges.
    gauges_dt : float
        Time interval for the gauges [s].

    Returns
    -------
    None
    """
    # create wave input
    spec_array = np.array([[1 / wave_period, wave_amplitude, 0]])

    # create the steering file
    create_steering_file(
        name_run,
        scheme,
        order_reconstruction,
        order_time_integration,
        courant_number,
        grid_size,
        domain_size,
        water_depth,
        run_time,
        output_interval,
        gauges_locations,
        gauges_dt,
        spec_array,
    )


    # Compile c++ code
    subprocess.run(["make"], cwd="code", stdout=subprocess.DEVNULL)

    # Run the model
    print(f"Running {name_run}...")
    if subprocess.run(["./barracuda", f"inputs/{name_run}.steer"], cwd="code").returncode != 0:
        return


    # Read the results
    print(f"\n")
    print(f"Reading {name_run} results and converting to xarray.Dataset...")
    ds = read_results(name_run)

    ds.attrs['wave_amplitude'] = wave_amplitude
    ds.attrs['wave_period'] = wave_period

    # Delete the results from the code/results folder
    subprocess.run(["rm", "-r", f"code/results/{name_run}"], cwd=".")
    subprocess.run(["rm", f"code/inputs/{name_run}.steer"], cwd=".")

    print(f"Done!")

    return ds



def run_spectral_wave(
    name_run: str,
    scheme: scheme,
    order_reconstruction: int,
    order_time_integration: int,
    courant_number: float,
    grid_size: float,
    domain_size: float,
    water_depth: float,
    spectrum_type: spec,
    significant_wave_height: float,
    peak_wave_period: float,
    run_time: float,
    output_interval: float,
    gauges_locations: np.ndarray,
    gauges_dt: float,
) -> None:
    """Run the wave model with a spectral wave input.

    Parameters
    ----------
    name_run : str
        Name of the run.
    scheme : scheme
        Numerical scheme used in the model.
    order_reconstruction : int
        Order of the flux reconstruction.
    order_time_integration : int
        Order of the time integration.
    courant_number : float
        Courant number.
    grid_size : float
        Grid size [m].
    domain_size : float
        Domain size [m].
    water_depth : float
        Water depth [m].
    spectrum_type : spec
        Type of empirical spectrum.
    significant_wave_height : float
        Significant wave height [m].
    peak_wave_period : float
        Peak period [s].
    run_time : float
        Simulation run time [s].
    output_interval : float
        Output interval [s].
    gauges_locations : array_like [Ng]
        Locations of the gauges [m]. 
        The reference point is the center of the wavemaker.
        Ng is the number of gauges.
    gauges_dt : float
        Time interval for the gauges [s].

    Returns
    -------
    None
    """

    # create wave input
    spec_all = create_empirical_spec(
        spectrum_type, significant_wave_height, peak_wave_period, water_depth, run_time
    )
    spec_array = spec_all[:, [0, 2, 3]]

    # create the steering file
    create_steering_file(
        name_run,
        scheme,
        order_reconstruction,
        order_time_integration,
        courant_number,
        grid_size,
        domain_size,
        water_depth,
        run_time,
        output_interval,
        gauges_locations,
        gauges_dt,
        spec_array,
    )

    # Compile c++ code
    subprocess.run(["make"], cwd="code", stdout=subprocess.DEVNULL)

    # Run the model
    print(f"Running {name_run}...")

    if subprocess.run(["./barracuda", f"inputs/{name_run}.steer"], cwd="code").returncode != 0:
        return

    # Read the results
    print(f"\n")
    print(f"Reading {name_run} results and converting to xarray.Dataset...")
    ds = read_results(name_run)

    ds.attrs['significant_wave_height'] = significant_wave_height
    ds.attrs['peak_wave_period'] = peak_wave_period


    # Delete the results from the code/results folder
    subprocess.run(["rm", "-r", f"code/results/{name_run}"], cwd=".")
    subprocess.run(["rm", f"code/inputs/{name_run}.steer"], cwd=".")

    print(f"Done!")

    return ds


def create_steering_file(
    name_run: str,
    scheme: scheme,
    order_reconstruction: int,
    order_time_integration: int,
    courant_number: float,
    grid_size: float,
    domain_size: float,
    water_depth: float,
    run_time: float,
    output_interval: float,
    gauges_locations: np.ndarray,
    gauges_dt: float,
    spec_array: np.ndarray,
) -> None:
    """Create the steering file for the model run and save it in code/inputs.

    Parameters
    ----------
    name_run : str
        Name of the run.
    scheme : scheme
        Numerical scheme used in the model.
    order_reconstruction : int
        Order of the flux reconstruction.
    order_time_integration : int
        Order of the time integration.
    courant_number : float
        Courant number.
    grid_size : float
        Grid size [m].
    domain_size : float
        Domain size [m].
    water_depth : float
        Water depth [m].
    run_time : float
        Simulation run time [s].
    output_interval : float
        Output interval [s].
    gauges_locations : array_like [Ng]
        Locations of the gauges [m]. 
        The reference point is the center of the wavemaker.
        Ng is the number of gauges.
    gauges_dt : float
        Time interval for the gauges [s].
    spec_array : array_like [Nf, 3]

    Returns
    -------
    None
    """
    # directory to save the steering file
    if not os.path.exists("code/inputs"):
        os.makedirs("code/inputs")

    gauges_locations = np.array(gauges_locations)
    for i in range(len(gauges_locations)):
        gauges_locations[i] = float(gauges_locations[i])
    gauges_locations_str = np.array2string(gauges_locations, separator=', ', max_line_width=np.inf)

    with open(f"code/inputs/{name_run}.steer", "w") as f:
        f.write(f"SCHEME: {scheme.name}\n")
        f.write(f"ORDER_RECONSTRUCTION: {order_reconstruction}\n")
        f.write(f"ORDER_TIME_INTEGRATION: {order_time_integration}\n")
        f.write(f"COURANT_NUMBER: {courant_number}\n")
        f.write(f"GRID_SIZE: {grid_size}\n")
        f.write(f"DOMAIN_SIZE: {domain_size}\n")
        f.write(f"WATER_DEPTH: {water_depth}\n")
        f.write(f"RUN_TIME: {run_time}\n")
        f.write(f"OUTPUT_INTERVAL: {output_interval}\n")
        f.write(f"GAUGES_LOCATIONS: {gauges_locations_str}\n")
        f.write(f"GAUGES_DT: {gauges_dt}\n")
        f.write("WAVE_INPUT_TABLE: FREQUENCY AMPLITUDE PHASE\n")
        for row in spec_array:
            f.write(f"{row[0]} {row[1]} {row[2]}\n")

        f.flush()
        os.fsync(f.fileno())


def read_output(name_run: str) -> np.ndarray:
    return np.zeros((1, 1))




def run_suzuki(
    name_run: str,
    scheme: scheme,
    order_reconstruction: int,
    order_time_integration: int,
    courant_number: float,
    grid_size: float,
    domain_size: float,
    water_depth: float,
    spec_all : np.ndarray,
    run_time: float,
    output_interval: float,
    gauges_locations: np.ndarray,
    gauges_dt: float,
) -> None:
    """Run the wave model with a spectral wave input.

    Parameters
    ----------
    name_run : str
        Name of the run.
    scheme : scheme
        Numerical scheme used in the model.
    order_reconstruction : int
        Order of the flux reconstruction.
    order_time_integration : int
        Order of the time integration.
    courant_number : float
        Courant number.
    grid_size : float
        Grid size [m].
    domain_size : float
        Domain size [m].
    water_depth : float
        Water depth [m].
    spectrum_type : spec
        Type of empirical spectrum.
    significant_wave_height : float
        Significant wave height [m].
    peak_wave_period : float
        Peak period [s].
    run_time : float
        Simulation run time [s].
    output_interval : float
        Output interval [s].
    gauges_locations : array_like [Ng]
        Locations of the gauges [m]. 
        The reference point is the center of the wavemaker.
        Ng is the number of gauges.
    gauges_dt : float
        Time interval for the gauges [s].

    Returns
    -------
    None
    """

    # create wave input
    # spec_all = create_empirical_spec(
    #     spectrum_type, significant_wave_height, peak_wave_period, water_depth, run_time
    # )
    spec_array = spec_all[:, [0, 2, 3]]

    # create the steering file
    create_steering_file(
        name_run,
        scheme,
        order_reconstruction,
        order_time_integration,
        courant_number,
        grid_size,
        domain_size,
        water_depth,
        run_time,
        output_interval,
        gauges_locations,
        gauges_dt,
        spec_array,
    )

    # Compile c++ code
    subprocess.run(["make"], cwd="code", stdout=subprocess.DEVNULL)

    # Run the model
    print(f"Running {name_run}...")

    if subprocess.run(["./barracuda", f"inputs/{name_run}.steer"], cwd="code").returncode != 0:
        return

    # Read the results
    print(f"\n")
    print(f"Reading {name_run} results and converting to xarray.Dataset...")
    ds = read_results(name_run)

    # ds.attrs['significant_wave_height'] = significant_wave_height
    # ds.attrs['peak_wave_period'] = peak_wave_period


    # Delete the results from the code/results folder
    subprocess.run(["rm", "-r", f"code/results/{name_run}"], cwd=".")
    subprocess.run(["rm", f"code/inputs/{name_run}.steer"], cwd=".")

    print(f"Done!")

    return ds