"""
Reorganize the checkboxes, such that subtype is underneath type in tree format.
"""
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dte

from .callbacks import available_checkboxes
from .data import dropdown_headers, dropdown_columns
from .components import Col, Link, Row


dr_types = ('UMAP', 't-SNE')

graph_div = html.Div([

        # container for 3D scatter plot
        dcc.Graph(
            id='plots-3d-scatter',
            figure={},
        ),

        # dropdown menu for umap/tsne.
        # wrapped in a row so that I could right align it.
        Row([
            html.Div([], className='col-auto mr-auto'),
            html.Div([
                dcc.Dropdown(
                    id='dr-type',
                    options=[{'label': label, 'value': label.lower()}
                              for label in dr_types],
                    value=dr_types[0].lower(),
                    clearable=False,
                    searchable=False,
                )], style={'width': '120px'}, className='col-auto',
            ),
        ]),
    ], className='col-xs-12 col-sm-7', id='left', #size='8',
)

control_divs = html.Div([

        Row([
            # Checkboxes for tumor or normal
            html.Div([
                html.H4("Show:"),
            ], className='col-sm-6'),
            html.Div([
                dcc.Checklist(
                    id='tumor-or-normal-checkboxes',
                    options=[{'label': 'Tumor', 'value': 'tumor'},
                             {'label': 'Matched Normal', 'value': 'normal'}],
                    values=['tumor', 'normal'],
                    inputStyle={'margin-right': '5px', 'margin-left': '-20px'},
                    labelStyle={'display': 'block',
                                'margin-bottom': '0',
                                'margin-left': '20px', 'font-size': '0.9rem',
                                'white-space': 'pre-wrap',
                            },
                )], className='col-sm-6'),
            ],
        ),

        html.Hr(),
        Row([
            # Radio for coloring
            html.Div([
                html.H4("Color samples by:"),
            ], className='col-sm-6'),
            html.Div([
                dcc.RadioItems(
                    id='color-radio',
                    options=[{'label': 'Tumor type', 'value': 'tumor.type'},
                             {'label': 'Organ system', 'value': 'organ_system'},
                             {'label': 'Subtype', 'value': 'subtype_full'}],
                    value='tumor.type',
                    inputStyle={'margin-right': '5px', 'margin-left': '-20px'},
                    labelStyle={'display': 'block',
                                'margin-bottom': '0',
                                'margin-left': '20px', 'font-size': '0.9rem',
                                'white-space': 'pre-wrap',
                            },
                )], className='col-sm-6'),
            ],
        ),

        html.Hr(),
        html.Div([
            # Multi Dropdown for each
            html.Div([
                html.H4("Specific tumor types:"),
            ]),
            dcc.Dropdown(
                id='dropdown-tumor-selection',
                options=[dict(disabled=False, **d) for d in \
                         available_checkboxes[dropdown_columns[0]]],
                value=[d['value'] for d in \
                         available_checkboxes[dropdown_columns[0]]],
                multi=True,
                clearable=True,
                searchable=True,
            ),
        ]),

        html.Div([
            html.Div([
                html.H4("Specific organ types:"),
            ]),
            dcc.Dropdown(
                id='dropdown-organ-selection',
                options=[dict(disabled=False, **d) for d in \
                         available_checkboxes[dropdown_columns[1]]],
                value=[d['value'] for d in \
                         available_checkboxes[dropdown_columns[1]]],
                multi=True,
                clearable=True,
                searchable=True,
            ),
        ]),

        html.Div([
            html.Div([
                html.H4("Specific subtypes:"),
            ]),
            dcc.Dropdown(
                id='dropdown-subtype-selection',
                options=[dict(disabled=False, **d) for d in \
                         available_checkboxes[dropdown_columns[2]]],
                value=[d['value'] for d in \
                         available_checkboxes[dropdown_columns[2]]],
                multi=True,
                clearable=True,
                searchable=True,
            )
        ]),
    ],
    className="hidden-xs col-sm-5", id='right',
    style={'display': 'block', 'min-width': '460px'}
)


page_dr = Row([
    # graphs
    graph_div,

    # menus
    control_divs,
])
