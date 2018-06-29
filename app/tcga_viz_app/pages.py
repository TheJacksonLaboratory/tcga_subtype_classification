import dash_core_components as dcc
import dash_html_components as html

from .page_dr import page_dr
from .page_about import page_about
from .page_hm import page_hm
from .components import Col, Row


def page_not_found(pathname):
    return html.P("No page '{}'".format(pathname))
