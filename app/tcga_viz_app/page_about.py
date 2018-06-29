import os
import dash
import dash_core_components as dcc
import dash_html_components as html

from .server import server
from .data import read_data
from .callbacks import *

static_page = os.path.join(server.root_path,
                           server.config['STATIC_FOLDER'],
                           'about_page.md')
with open(static_page) as fin:
    about = fin.read()

page_about = html.Div([
    dcc.Markdown(children=about)
])



