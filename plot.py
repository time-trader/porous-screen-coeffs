import plotly.graph_objects as go


def plot_forces(title, forces_dict):
    fig = go.Figure()
    for key in forces_dict.keys():
        fig.add_trace(go.Scatter(x=forces_dict[key][0], y=forces_dict[key][1],
                                 mode='lines', name=key))
    fig.update_layout(
        title={
            'text': title,
            'y':0.95,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        xaxis_title="AoA [deg]",
        yaxis_title="Forces [N]")
    fig.show()
