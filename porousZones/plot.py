import os.path

import plotly.graph_objects as go
import pandas as pd


def plot_forces(title, forces_dict, save_eps=False):
    fig = go.Figure()


    for key in forces_dict.keys():
        if "Analytical" in key:
            _mode = 'lines'
        else:
            _mode = 'markers'

        fig.add_trace(go.Scatter(x=forces_dict[key][0], y=forces_dict[key][1],
                                 mode=_mode, name=key,
                                 marker=dict(
                                    # size=5,
                                 ))
                      )

    fig.update_layout(
        title={
            #'text': title,
            'y':0.95,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        # xaxis=dict(range=[-90, 90]),
        xaxis_title="AoA [deg]",
        yaxis_title="Forces [N]",
        legend=dict(
            orientation="h",
            yanchor="top",
            y=1.25,
            xanchor="left",
            x=0.05),
        template='plotly_white',
        font_family="Computer Modern"
    )

    # fig.update_xaxes(dtick=30)

    if save_eps:
        fig.write_image(os.path.join("plots", "fig_{}.eps".format(title)), width=600, height=400, scale=20)
    else:
        fig.show()
