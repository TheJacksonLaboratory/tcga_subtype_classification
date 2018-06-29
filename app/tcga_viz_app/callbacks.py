import numpy as np
import time
from dash.dependencies import Output, Input, State, Event

from .server import app
from .data import *
from .plotting import generate_dr_figure, marginal_heatmap


sample_data = read_data()
available_checkboxes = get_available_checkboxes(sample_data)

# Put your callbacks here

@app.callback(
    Output('plots-3d-scatter', 'figure'),
    [
     Input('color-radio', 'value'),
     Input('tumor-or-normal-checkboxes', 'values'),
     Input('dropdown-tumor-selection', 'value'),
     Input('dropdown-organ-selection', 'value'),
     Input('dropdown-subtype-selection', 'value'),
     Input('dr-type', 'value')],
    [State('plots-3d-scatter', 'figure')]
    )
def update_dr_graph(color_value, tumor_or_normal, dropdown_tumor,
                    dropdown_organ, dropdown_subtype, dr_type,
                    figure):
    cols = ('tumor.type', 'organ_system', 'subtype_full')
    dropdowns = (dropdown_tumor, dropdown_organ, dropdown_subtype)

    subset = sample_data[sample_data.match.isin(tumor_or_normal)]
    for col, dropdown in zip(cols, dropdowns):
        if dropdown == []: continue
        subset = subset[subset[col].isin(dropdown)]

    subset_grouped = subset.groupby(('match',) + cols)

    layout = figure.get('layout', None)
    return generate_dr_figure(subset_grouped, color_value, layout, dr_type)

@app.callback(
    Output('plots-heatmap', 'figure'),
    [Input('dropdown-heatmap-choice', 'value')],
    [State('plots-heatmap', 'figure')
     ])
def update_primary_heatmap(data_choice, fig_state):
    if data_choice == 'primary_ext':
        return marginal_heatmap(*read_primary_ext_heatmap_data())
    return marginal_heatmap(*read_primary_heatmap_data(), margins=None)

def plot_subtype_heatmap(subtype, **kwargs):
    return marginal_heatmap(*read_subtype_heatmap_data(subtype), **kwargs)
