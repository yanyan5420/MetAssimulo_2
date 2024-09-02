import plotly.graph_objects as go
import numpy as np


def plot_all_metabolites(mixture_list, processed_data_dict, ppm_scale):
    fig = go.Figure()
    config = dict({'scrollZoom': True})

    for meta_name in mixture_list:
        fig.add_trace(go.Scatter(
            x=ppm_scale,
            y=processed_data_dict[meta_name],
            name=meta_name,
        ))

    fig.update_xaxes(title_text="ppm")
    fig.update_yaxes(title_text="intensity")
    fig['layout']['xaxis']['autorange'] = "reversed"
    fig.update_layout(
        title_text="Metabolite spectrum in mixture",
        template="plotly_white",
        height=500,
        legend=dict(
            # orientation='h',
            # yanchor='top',
            # y=-0.3,
            # xanchor='center',
            # x=0.5,
            itemclick="toggleothers",
            itemdoubleclick="toggle"
        )
    )
    return fig


def plot_mean_spectra(final_data_dict, ppm_scale):
    fig = go.Figure()
    config = dict({'scrollZoom': True})

    fig.add_trace(
        go.Scatter(
            x=ppm_scale,
            y=final_data_dict["replicate_mean"],
            mode="lines",
            line=dict(color='skyblue'),
        )
    )
    fig.update_xaxes(title_text="ppm")
    fig.update_yaxes(title_text="intensity")
    fig['layout']['xaxis']['autorange'] = "reversed"
    fig.update_layout(
        title_text="Mean of all replicates",
        template="plotly_white",
    )
    return fig


def plot_replicate_spectra(final_data_dict, select_replicate, ppm_scale):
    fig = go.Figure()
    config = dict({'scrollZoom': True})

    fig.add_trace(
        go.Scatter(
            x=ppm_scale,
            y=final_data_dict[select_replicate],
            mode="lines",
            line=dict(color='darkseagreen'),
        )
    )
    fig.update_xaxes(title_text="ppm")
    fig.update_yaxes(title_text="intensity")
    fig['layout']['xaxis']['autorange'] = "reversed"
    fig.update_layout(
        title_text="Simulated mixture spectra for replicates",
        template="plotly_white",
    )
    return fig


def plot_replicate_spectra_with_peak_shift(final_data_dict, select_replicate, ppm_scale):
    fig = go.Figure()
    config = dict({'scrollZoom': True})

    fig.add_trace(
        go.Scatter(
            x=ppm_scale,
            y=final_data_dict[select_replicate][1],
            mode="lines",
            line=dict(color='darkseagreen'),
        )
    )
    fig.update_xaxes(title_text="ppm")
    fig.update_yaxes(title_text="intensity")
    fig['layout']['xaxis']['autorange'] = "reversed"
    fig.update_layout(
        title_text="Simulated mixture spectra for replicates",
        template="plotly_white",
    )
    return fig


def plot_stacked_spectra(final_data_dict, ppm_scale, num_replicates, v_space):
    fig = go.Figure()
    config = dict({'scrollZoom': True})
    v = 0
    for i in range(num_replicates):
        temp_name = 'replicate_' + str(i + 1)
        y_intensity = np.array(final_data_dict[temp_name])
        fig.add_trace(go.Scatter(
            x=ppm_scale,
            y=y_intensity + v,
            name=temp_name,
            hovertext=tuple(zip(ppm_scale, y_intensity)),
            hoverinfo="text",
        ))
        v = v + v_space

    fig.update_xaxes(title_text="ppm")
    fig.update_yaxes(title_text="intensity")
    fig.update_yaxes(tickvals=[])

    fig['layout']['xaxis']['autorange'] = "reversed"
    fig.update_layout(
        title_text="All replicates in the group",
        template="plotly_white",
        height=500,
        legend=dict(
            orientation='h',
            yanchor='top',
            y=-0.15,
            xanchor='center',
            x=0.5,
            itemclick="toggleothers",
            itemdoubleclick="toggle"
        )
    )
    return fig


def plot_stacked_spectra_with_same_ph(final_data_dict, ppm_scale, num_replicates, v_space):
    fig = go.Figure()
    config = dict({'scrollZoom': True})
    v = 0
    for i in range(num_replicates):
        temp_name = 'replicate_' + str(i + 1)
        ph, y = final_data_dict[temp_name]
        y_intensity = np.array(y)
        fig.add_trace(go.Scatter(
            x=ppm_scale,
            y=y_intensity + v,
            name=temp_name,
            hovertext=tuple(zip(ppm_scale, y_intensity)),
            hoverinfo="text",
        ))
        v = v + v_space

    fig.update_xaxes(title_text="ppm")
    fig.update_yaxes(title_text="intensity")
    fig.update_yaxes(tickvals=[])

    fig['layout']['xaxis']['autorange'] = "reversed"
    fig.update_layout(
        title_text="All replicates in the group",
        template="plotly_white",
        height=500,
        legend=dict(
            orientation='h',
            yanchor='top',
            y=-0.15,
            xanchor='center',
            x=0.5,
            itemclick="toggleothers",
            itemdoubleclick="toggle"
        )
    )
    return fig


def plot_stacked_spectra_with_diff_ph(final_data_dict, ppm_scale, num_replicates, v_space):
    fig = go.Figure()
    config = dict({'scrollZoom': True})
    v = 0
    annotations = []
    for i in range(num_replicates):
        temp_name = 'replicate_' + str(i + 1)
        ph, y = final_data_dict[temp_name]
        temp_ph = float(ph)
        y_intensity = np.array(y)
        fig.add_trace(go.Scatter(
            x=ppm_scale,
            y=y_intensity + v,
            name=temp_name + ", " + "pH=" + str(temp_ph),
            hovertext=tuple(zip(ppm_scale, y_intensity)),
            hoverinfo="text",
        ))
        # add labels for each line
        annotations.append(dict(xref="paper", x=1.0, y=0 + v,
                                xanchor="left", yanchor="middle",
                                text="pH=" + str(temp_ph),
                                font=dict(family='Arial', size=14),
                                showarrow=False))
        v = v + v_space

    fig.update_xaxes(title_text="ppm")
    fig.update_yaxes(title_text="intensity")
    fig.update_yaxes(tickvals=[])

    fig['layout']['xaxis']['autorange'] = "reversed"
    fig.update_layout(
        annotations=annotations,
        title_text="All replicates in the group",
        template="plotly_white",
        height=500,
        legend=dict(
            orientation='h',
            yanchor='top',
            y=-0.15,
            xanchor='center',
            x=0.5,
            itemclick="toggleothers",
            itemdoubleclick="toggle",
            traceorder="reversed"
        )
    )
    return fig


def plot_stacked_spectra_with_ph(sort_data_dict, ppm_scale, v_space):
    fig = go.Figure()
    config = dict({'scrollZoom': True})
    v = 0
    annotations = []
    for repli_name, [ph, y] in sort_data_dict.items():
        temp_ph = float(ph)
        y_intensity = np.array(y)
        fig.add_trace(go.Scatter(
            x=ppm_scale,
            y=y_intensity + v,
            name=repli_name + ", " + "pH=" + str(temp_ph),
            hovertext=tuple(zip(ppm_scale, y_intensity)),
            hoverinfo="text",
        ))
        # add labels for each line
        annotations.append(dict(xref="paper", x=1.0, y=0 + v,
                                xanchor="left", yanchor="middle",
                                text="pH=" + str(temp_ph),
                                font=dict(family='Arial', size=14),
                                showarrow=False))
        v = v + v_space

    fig.update_xaxes(title_text="ppm")
    fig.update_yaxes(title_text="intensity")
    fig.update_yaxes(tickvals=[])

    fig['layout']['xaxis']['autorange'] = "reversed"
    fig.update_layout(
        annotations=annotations,
        title_text="All replicates in the group",
        template="plotly_white",
        height=500,
        legend=dict(
            orientation='h',
            yanchor='top',
            y=-0.15,
            xanchor='center',
            x=0.5,
            itemclick="toggleothers",
            itemdoubleclick="toggle",
            traceorder="reversed"
        )
    )
    return fig
