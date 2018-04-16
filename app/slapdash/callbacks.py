import numpy as np
import time
from dash.dependencies import Output, Input, State, Event

from .server import app
from .data import *
from .plotting import generate_figure



sample_data = read_data()
available_checkboxes = get_available_checkboxes(sample_data)


# Put your callbacks here


@app.callback(
    Output('tumortype-checkboxes', 'options'),
    [Input('dropdown-data-choice', 'value')])
def set_checkbox_options1(header_value):
    disabled = header_value != 'tumor.type.matched'
    options = [{'disabled': disabled, **val}
              for val in available_checkboxes['tumor.type.matched']]
    return options

@app.callback(
    Output('cancertype-checkboxes', 'options'),
    [Input('dropdown-data-choice', 'value')])
def set_checkbox_options2(header_value):
    disabled = header_value != 'cancertype'
    options = [{'disabled': disabled, **val}
              for val in available_checkboxes['cancertype']]
    return options

@app.callback(
    Output('cancersubtype-checkboxes', 'options'),
    [Input('dropdown-data-choice', 'value')])
def set_checkbox_options3(header_value):
    disabled = header_value != 'subtype'
    options = [{'disabled': disabled, **val}
              for val in available_checkboxes['subtype']]
    return options


@app.callback(
    Output('tumortype-checkboxes', 'values'),
    [Input('tumortype-checkboxes', 'options'),
     Input('clear-button-count', 'children')])
def set_checkbox_values1(options, count_children):
    print(count_children)
    counter = int(count_children[0])
    if counter > 0:
        return []
    return [opt['value'] for opt in options if not opt['disabled']]

@app.callback(
    Output('cancertype-checkboxes', 'values'),
    [Input('cancertype-checkboxes', 'options')])
def set_checkbox_values2(options):
    return [opt['value'] for opt in options if not opt['disabled']]

@app.callback(
    Output('cancersubtype-checkboxes', 'values'),
    [Input('cancersubtype-checkboxes', 'options')])
def set_checkbox_values3(options):
    return [opt['value'] for opt in options if not opt['disabled']]


@app.callback(
    Output('clear-button-count', 'children'),
    [Input('clear-tumor-button', 'n_clicks'),
     Input('dropdown-data-choice', 'value')])
def update_counter(n_clicks, value):
    print('updating counter', n_clicks, value)
    if n_clicks is None or value != 'tumor.type.matched':
        print('this path')
        return '0'
    return '1'



@app.callback(
    Output('plots-3d-scatter', 'figure'),
    [Input('tumortype-checkboxes', 'options'),
     Input('cancertype-checkboxes', 'options'),
     Input('cancersubtype-checkboxes', 'options'),
     Input('tumortype-checkboxes', 'values'),
     Input('cancertype-checkboxes', 'values'),
     Input('cancersubtype-checkboxes', 'values')
     ])
def update_graph(tumor_options, type_options, subtype_options,
                tumor_values, type_values, subtype_values):
    options = (tumor_options, type_options, subtype_options)
    values = (tumor_values, type_values, subtype_values)
    cols = ('tumor.type.matched', 'cancertype', 'subtype')
    ind = np.argmax(list(map(len, values)))

    return generate_figure(sample_data, options[ind], values[ind], cols[ind])
