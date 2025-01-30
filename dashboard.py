from dash import Dash, dash_table, dcc, html
from dash.dependencies import Input, Output, State
import sys, os
from pymatgen.analysis.wulff import WulffShape
from pymatgen.core.structure import Lattice
from pymatgen.ext.matproj import MPRester
mpr = MPRester('mlcC4gtXFVqN9WLv')

from scipy.spatial.qhull import QhullError

app = Dash(__name__)
server = app.server


app.layout = html.Div([
    dcc.Tabs([
        dcc.Tab([
            # make a table for the abc lattice parameters
            dash_table.DataTable(
                id='abc',
                columns=[{'name': 'a', 'id': 'a',
                          'deletable': False, 'renamable': False},
                         {'name': 'b', 'id': 'b',
                          'deletable': False, 'renamable': False},
                         {'name': 'c', 'id': 'c',
                          'deletable': False, 'renamable': False},
                         ],
                data=[{'a': 1, 'b': 1, 'c': 1}],
                editable=True,
                row_deletable=False,
                style_cell={"textAlign": "center", 'minWidth': '100px'}, 
                fill_width=False
            ),

            # make a table for the angles lattice parameters
            dash_table.DataTable(
                id='angles',
                columns=[{'name': 'alpha', 'id': 'alpha',
                          'deletable': False, 'renamable': False},
                         {'name': 'beta', 'id': 'beta',
                          'deletable': False, 'renamable': False},
                         {'name': 'gamma', 'id': 'gamma',
                          'deletable': False, 'renamable': False},
                         ],
                data=[{'alpha': 90, 'beta': 90, 'gamma': 90}],
                editable=True,
                row_deletable=False,
                style_cell={ "textAlign": "center", 'minWidth': '100px'}, 
                fill_width=False
            ),

            # make a table for the miller index facets and surface energy
            dash_table.DataTable(
                id='hkl_and_surface_energy',
                columns=[{'name': 'h', 'id': 'h', 'deletable': False, 'renamable': False},
                         {'name': 'k', 'id': 'k', 'deletable': False, 'renamable': False},
                         {'name': 'l', 'id': 'l', 'deletable': False, 'renamable': False},
                         {'name': 'Surface energy', 'id': 'surface_energy', 
                          'deletable': False, 'renamable': False},
                         {'name': 'Area fraction', 'id': 'area_frac', 
                          'deletable': False, 'renamable': False}],
                data=[{'h': 1, 'k': 0, 'l': 0, 'surface_energy': 1, 'area_frac': 1}],
                editable=True,
                row_deletable=True,
                style_cell={"textAlign": "center", 'minWidth': '100px'},
                fill_width=False
            ),

            # add a button for adding more facets
            html.Button('Add surface', id='editing-rows-button', n_clicks=0),
            # add a box for inputting specific mpid
            dcc.Input(id="MPID", type="text", placeholder="MPID", style={'marginRight':'10px'}, debounce=True),
            dcc.Graph(id='wulff_shape'),

            dcc.Upload(id='slab_vrun',
                       children=html.Div(['Slab ', html.A('vasprun.xml file')]),
                       style={'width': '100%', 'height': '60px', 'lineHeight': '60px',
                              'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                              'textAlign': 'center', 'margin': '10px'}),
            dcc.Upload(id='bulk_vrun',
                       children=html.Div(['Bulk ', html.A('vasprun.xml file')]),
                       style={'width': '100%', 'height': '60px', 'lineHeight': '60px',
                              'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                              'textAlign': 'center', 'margin': '10px'}),
            # make a table for the miller index facets for uploaded xml files
            dash_table.DataTable(id='hkl_xml', columns=[{'name': 'h', 'id': 'h', 'deletable': False, 'renamable': False},
                                                        {'name': 'k', 'id': 'k', 'deletable': False, 'renamable': False},
                                                        {'name': 'l', 'id': 'l', 'deletable': False, 'renamable': False}],
                                 data=[{'h': 1, 'k': 0, 'l': 0}], editable=True, row_deletable=False, 
                                 style_cell={"textAlign": "center", 'minWidth': '100px'}, fill_width=False),
            html.Button('Calculate surface energy', id='calculate_button', n_clicks=0),
        ], label='Wulff Shapes')
    ])
])

if __name__ == '__main__':
    app.run_server(debug=True)
