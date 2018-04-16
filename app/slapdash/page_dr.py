"""
Reorganize the checkboxes, such that everything can be selected at once.
"""
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dte

from .callbacks import *
from .components import Col, Link, Row


graph_div = html.Div([
        dcc.Graph(
            id='plots-3d-scatter',
            figure={},
        )
    ], className='col-xs-12 col-sm-6', id='left', #size='8',
)

control_divs = Col([
        html.Div([
            dcc.Dropdown(
                id='dropdown-data-choice',
                options=[{'label': header, 'value': col} for header, col in \
                         zip(dropdown_headers, dropdown_columns)],
                value=dropdown_columns[0],
            )],
        ),

        html.Div([
            dcc.Checklist(
                id='dropdown-checkboxes',
                options=[{'label': val, 'value': val}
                         for val in available_checkboxes[dropdown_columns[0]]],
                values=available_checkboxes[dropdown_columns[0]],
                inputStyle=dict(margin='5'),
                labelStyle={'display': 'block', 'margin-bottom': '0'},
            )], style=dict(display='block', columnCount=2)
        ),

    ],
)


control_divs = html.Div([
        html.Div([
            html.Div([
                html.H3("Select filter property:"),
            ]),

            dcc.Dropdown(
                id='dropdown-data-choice',
                options=[{'label': header, 'value': col} for header, col in \
                         zip(dropdown_headers, dropdown_columns)],
                value=dropdown_columns[0],
            )],
        ),

        html.Hr(),
        html.Div([
            html.H4("Tumor Type"),
            html.Button('Clear selection', id='clear-tumor-button'),
            html.Div(id='clear-button-count', style={'display': 'none'}
            ),
        ]),

        html.Div([
            dcc.Checklist(
                id='tumortype-checkboxes',
                options=[dict(disabled=False, **d) for d in \
                         available_checkboxes[dropdown_columns[0]]],
                values=[d['value'] for d in \
                        available_checkboxes[dropdown_columns[0]]],
                inputStyle={'margin-right': '5px', 'margin-left': '-20px'},
                labelStyle={'display': 'block', 'margin-bottom': '0', 'margin-left': '20px'},
            )], style=dict(display='block', columnCount=4)
        ),

        html.Hr(),
        html.Div([
            html.H4("Cancer Type"),
        ]),

        html.Div([
            dcc.Checklist(
                id='cancertype-checkboxes',
                options=[dict(disabled=True, **d) for d in \
                         available_checkboxes[dropdown_columns[1]]],
                values=[d['value'] for d in \
                        available_checkboxes[dropdown_columns[1]]],
                inputStyle={'margin-right': '5px', 'margin-left': '-20px'},
                labelStyle={'display': 'block', 'margin-bottom': '0', 'margin-left': '20px'},
            )], style=dict(display='block', columnCount=2)
        ),

        html.Hr(),
        html.Div([
            html.H4("Cancer SubType"),
        ]),

        html.Div([
            dcc.Checklist(
                id='cancersubtype-checkboxes',
                options=[dict(disabled=True, **d) for d in \
                         available_checkboxes[dropdown_columns[2]]],
                values=[d['value'] for d in \
                        available_checkboxes[dropdown_columns[2]]],
                inputStyle={'margin-right': '5px', 'margin-left': '-20px'},
                labelStyle={'display': 'block', 'margin-bottom': '0', 'margin-left': '20px'},
            )], style=dict(display='block', columnCount=2)
        ),
        ],
        className="hidden-xs col-sm-6", id='right',
        style={'display': 'block'}
)


page_dr = Row([
    # graphs
    graph_div,

    # menus
    control_divs,
])
