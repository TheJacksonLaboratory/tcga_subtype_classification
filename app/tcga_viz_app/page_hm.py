import dash
import dash_core_components as dcc
import dash_html_components as html
#import dash_table_experiments as dte

import numpy as np
import pandas as pd

from .callbacks import available_checkboxes, plot_subtype_heatmap
from .data import dropdown_headers, dropdown_columns
from .components import Col, Link, Row


container = html.Div([
    html.Div([
        html.H4('Primary classification results'),
    ]),

    html.Div([
        dcc.Dropdown(
            id='dropdown-heatmap-choice',
            options=[{'label': lab, 'value': val} for lab, val in \
                     zip(('Primary Cross-validation', 'Primary External validation'),
                         ('primary_cv', 'primary_ext'))],
            value='primary_cv',
            clearable=False,
            searchable=False,
        )],
    ),

    html.Div([
        dcc.Graph(
            id='plots-heatmap'
        ),
    ]),

    html.Div([
        html.H4('Subtype classification results'),
    ]),

    html.Div([], id='subtype-plot-container')
], className="container")

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]
subtypes = ['BRCA', 'HNSC', 'KIRC', 'KIRP', 'LGG', 'LUAD', 'LUSC', 'OV', 'PRAD',
        'SKCM', 'STAD']

def subtype_graph(subtype, size=6, offset=0):
    graph = dcc.Graph(id=f'plots-heatmap-{subtype}')
    graph.figure = plot_subtype_heatmap(subtype, size=None,
                                        margins=dict(t=150, b=10, l=200, r=10),
                                        title=None)#subtype)
    className = f"col-sm-{size}"
    if offset > 0:
        className += f" col-sm-offset-{offset}"
    return html.Div([html.Div([html.H5(subtype)], style={'text-align': 'center'}),
                     graph, ], className=className)

for subtype_subset in chunks(subtypes, 2):
    offset = 0
    if len(subtype_subset) < 2:
        offset = (12 - 6*len(subtype_subset)) // 2
        row = [subtype_graph(subtype_subset[0], offset=offset)]
               #subtype_graph(subtype_subset[1], offset=0)]
    else:
        row = [subtype_graph(subtype, offset=offset) for subtype in subtype_subset]
    container.children[-1].children.append(Row(row, style={'margin': '100px 25px'}))

page_hm = Row([
    container,
])


#table_cols = ['sampleid', 'patientid', 'sample.type', 'tumor.type', 'subtype',
#        'cancertype']
#table_headers = ['ID', 'Patient ID', 'Sample type', 'Tumor type', 'Subtype',
#        'Cancer Type']
#table_div = html.Div([
#    dte.DataTable(
#        rows=sample_data.loc[:,table_cols].to_dict('records'),
#        columns=table_cols,
#        row_selectable=True,
#        filterable=True,
#        sortable=True,
#        selected_row_indices=[0],
#        id='datatable',
#    )],
#    style=dict(width='100%', display='block')
#)
