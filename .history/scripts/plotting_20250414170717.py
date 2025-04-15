"""
This module contains functions to plot the numerical results.

It includes functions to create animated plots of the free surface evolution and to plot the free surface time series at specified gauges.
The module uses Plotly for interactive visualizations and provides options for plotting both numerical and exact solutions.

Functions
---------
create_animation_plotly :
    Plots the free surface evolution over time, creating an animated figure with optional exact solution data.

plot_gauges :
    Plots the free surface evolution at specified gauges, with the option to display both numerical results and the exact solution.

"""

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import xarray as xr
from plotly.subplots import make_subplots
import os


def create_animation_plotly(
    ds_list: list, ds_exact: xr.Dataset = None, frame_duration: int = 200
) -> go.Figure:
    """Plot the free surface evolution.

    Parameters
    ----------
    ds_list : array_like [xarray.Dataset]
        Dataset with the numerical results.

    ds_exact : xarray.Dataset, optional
        Dataset with the exact solution.

    frame_duration : int
        The duration between each frame in milliseconds.
        default: 200
        To speed up the animation, increase this value.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Figure with the free surface evolution.
    """

    ds0 = ds_list[0]
    time = ds0.time.values

    # create figure
    fig = go.Figure()

    # Add trace for exact solution
    if ds_exact is not None:
        fig.add_trace(
            go.Scatter(
                x=ds_exact.x.values,
                y=ds_exact.eta[0, :],
                mode="lines",
                line=dict(width=3, color="black"),
                name="Exact solution",
                opacity=0.4,
            ),
        )

    # Add traces for initial frame
    for ds in ds_list:
        fig.add_trace(
            go.Scatter(
                x=ds.x.values,
                y=ds.eta[0, :],
                mode="lines",
                line=dict(width=3),
                name=ds.attrs["name_run"],
            ),
        )

    # create frames
    if ds_exact is not None:
        frames = [
            dict(
                name=k,
                data=[
                    go.Scatter(x=ds_exact.x.values, y=ds_exact.eta[k, :], mode="lines")
                ]
                + [
                    go.Scatter(x=ds.x.values, y=ds.eta[k, :], mode="lines")
                    for ds in ds_list
                ],
                traces=[i for i in range(len(ds_list) + 1)],
            )
            for k in range(len(time))
        ]
    else:
        frames = [
            dict(
                name=k,
                data=[
                    go.Scatter(x=ds.x.values, y=ds.eta[k, :], mode="lines")
                    for ds in ds_list
                ],
                traces=[i for i in range(len(ds_list))],
            )
            for k in range(len(time))
        ]

    # create buttons
    updatemenus = [
        dict(
            type="buttons",
            buttons=[
                dict(
                    label="&#9658;",
                    method="animate",
                    args=[
                        [f"{k}" for k in range(len(time))],
                        dict(
                            mode="immediate",
                            frame=dict(duration=frame_duration, redraw=True),
                            transition=dict(duration=0),
                        ),
                    ],
                ),
                dict(
                    label="&#9209;",
                    method="animate",
                    args=[
                        [None],
                        dict(
                            mode="immediate",
                            frame=dict(duration=0, redraw=True),
                            transition=dict(duration=0),
                        ),
                    ],
                ),
                dict(
                    label="&#9198;",
                    method="animate",
                    args=[
                        [0],
                        dict(
                            mode="immediate",
                            frame=dict(duration=0, redraw=True),
                            transition=dict(duration=0),
                        ),
                    ],
                ),
            ],
            showactive=False,
        ),
    ]

    # create sliders
    sliders = [
        {
            "currentvalue": {
                "prefix": "<b>Time: ",
                "visible": True,
                "xanchor": "right",
                "font": {"size": 16},
            },
            "steps": [
                {
                    "args": [
                        [k],
                        dict(
                            mode="immediate",
                            frame=dict(duration=0, redraw=True),
                            transition=dict(duration=0),
                        ),
                    ],
                    "label": str(np.around(time[k], decimals=1)),
                    "method": "animate",
                }
                for k in range(len(time))
            ],
        }
    ]

    # update layout
    fig.update(frames=frames),
    fig.update_layout(
        updatemenus=updatemenus,
        sliders=sliders,
        template="plotly_white",
    )

    # add labels
    fig.update_layout(
        xaxis_title="x [m]",
        yaxis_title="Free surface elevation [m]",
    )

    fig.update_layout(font=dict(family="monospace"))

    return fig


def plot_gauges(
    ds_list: list, ds_exact: xr.Dataset = None, gauges: list = None
) -> go.Figure:
    """Plot the free surface evolution at the gauges.

    Parameters
    ----------
    ds_list : array_like [xarray.Dataset]
        Dataset with the numerical results.

    ds_exact : xarray.Dataset, optional
        Dataset with the exact solution.

    gauges : list, optional
        List with the gauges to plot. If None, all gauges are plotted.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Figure with the free surface evolution at the different gauges.
    """

    # color palette
    colors = ["blue", "red", "green", "cyan", "magenta", "yellow", "black", "orange"]

    if gauges is not None:
        index0 = []  # index of the gauge inside the gauges array
        index_x = []  # index of the gauge inside the x array
        for ds in ds_list + [ds_exact]:
            new_index = np.argmin(
                np.abs(ds.x.values[ds.index_gauges.values][:, np.newaxis] - gauges),
                axis=0,
            )
            index0.append(new_index)
            index_x.append(ds.index_gauges.values[new_index])

        # create subplot for each gauge
        subplot_title = [
            f"G{i+1} = {ds_list[0].x[index_x[0][i]]:.2f} m"
            for i in range(len(index0[0]))
        ]
        fig = make_subplots(
            rows=len(index0[0]), cols=1, shared_xaxes=True, subplot_titles=subplot_title
        )

        # loop over the gauges
        for num, i in enumerate(index0[0]):

            if ds_exact is not None:
                fig.add_trace(
                    go.Scatter(
                        x=ds_exact.time_gauges.values,
                        y=ds_exact.eta_gauges[:, i],
                        mode="lines",
                        line=dict(width=3, color="black"),
                        name="Exact solution",
                        opacity=0.4,
                    ),
                    row=num + 1,
                    col=1,
                )

            # Add traces for numerical results
            for j, ds in enumerate(ds_list):
                fig.add_trace(
                    go.Scatter(
                        x=ds.time_gauges.values,
                        y=ds.eta_gauges[:, i],
                        mode="lines",
                        line=dict(width=3, color=colors[j % len(colors)]),
                        name=ds.attrs["name_run"],
                    ),
                    row=num + 1,
                    col=1,
                )

            # add labels
            fig.update_yaxes(title_text="Free surface [m]", row=num + 1, col=1)
        fig.update_xaxes(title_text="Time [s]", row=num + 1, col=1)

        # update layout
        fig.update_layout(template="plotly_white")
        fig.update_layout(font=dict(family="monospace"))

        return fig

    subplot_title = [
        f"G{i+1} = {ds_list[0].x[ds_list[0].index_gauges[i]]:.2f} m"
        for i in range(ds_list[0].eta_gauges.shape[1])
    ]
    fig = make_subplots(
        rows=ds_list[0].eta_gauges.shape[1],
        cols=1,
        shared_xaxes=True,
        subplot_titles=subplot_title,
    )

    # loop over the gauges
    for i in range(ds_list[0].eta_gauges.shape[1]):

    for num, i in enumerate(index0[0]):
        
        # Add trace for exact solution
        if ds_exact is not None:
            fig.add_trace(
                go.Scatter(
                    x=ds_exact.time_gauges.values,
                    y=ds_exact.eta_gauges[:, i],
                    mode="lines",
                    line=dict(width=3, color="black"),
                    name="Exact solution",
                    opacity=0.4,
                ),
                row=i + 1,
                col=1,
            )

        # Add traces for numerical results
        for j, ds in enumerate(ds_list):
            fig.add_trace(
                go.Scatter(
                    x=ds.time_gauges.values,
                    y=ds.eta_gauges[:, i],
                    mode="lines",
                    line=dict(width=3, color=colors[j % len(colors)]),
                    name=ds.attrs["name_run"],
                ),
                row=i + 1,
                col=1,
            )

        # add labels
        fig.update_yaxes(title_text="Free surface [m]", row=num + 1, col=1)
    fig.update_xaxes(title_text="Time [s]", row=num + 1, col=1)

    # update layout
    fig.update_layout(template="plotly_white")
    fig.update_layout(font=dict(family="monospace"))

    # remove the xaxis from first N-1 gauges
    for i in range(ds_list[0].eta_gauges.shape[1] - 1):
        fig.update_xaxes(showticklabels=False, row=i + 1, col=1)
        fig.update_xaxes(title_text="", row=i + 1, col=1)

    return fig
