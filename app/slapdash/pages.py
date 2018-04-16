import dash_core_components as dcc
import dash_html_components as html

from .page_dr import page_dr
from .page_about import page_about
from .components import Col, Row


#page_dr = html.Div("Page 1")
page_hm = html.Div("Page 2")


def page_not_found(pathname):
    return html.P("No page '{}'".format(pathname))
