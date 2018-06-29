import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dte

from .callbacks import available_checkboxes
from .data import dropdown_headers, dropdown_columns
from .components import Col, Link, Row
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

# upload csv via form
# upload sample labels
# gene intersection
# model classification
# compare classification with provided labels
