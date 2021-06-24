import plotly.graph_objects as go
import pandas as pd

def plot_forces(title, forces_dict):
    fig = go.Figure()

    color_list = [
        "rgb(0, 0, 0)",
        "rgb(120, 120, 120)",
        "rgb(180, 180, 180)"
    ]

    marker_list = [3, 4, 5, 6]

    marker_i=0
    color_i=0

    for key in forces_dict.keys():
        if "Analytical" in key:
            line_dict = dict(color=color_list[color_i])
            _mode = 'lines'
            color_i+=1
        else:
            line_dict = dict()
            _mode = 'markers'
            marker_i+=1

        fig.add_trace(go.Scatter(x=forces_dict[key][0], y=forces_dict[key][1],
                                 mode=_mode, name=key,
                                 line=line_dict,
                                 marker=dict(
                                    symbol = marker_list[marker_i],
                                    size=15,
                                    color="black"
                                 ))
                      )

    fig.update_layout(
        title={
            'text': title,
            'y':0.95,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        xaxis_title="AoA [deg]",
        yaxis_title="Forces [N]",
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99),
        template='plotly_white'

    )
    fig.show()
