from dash import Dash, dash_table, dcc, html
from dash.dependencies import Input, Output, State
import sys
from pymatgen.analysis.wulff import WulffShape
from pymatgen.core.structure import Lattice
from pymatgen.ext.matproj import MPRester
mpr = MPRester('mlcC4gtXFVqN9WLv')

app = Dash(__name__)
server = app.server

app.layout = html.Div([

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
        style_cell={ "textAlign": "center", 'minWidth': '100px'}, 
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
        columns=[{'name': 'h', 'id': 'h',
                  'deletable': False, 'renamable': False},
                 {'name': 'k', 'id': 'k',
                  'deletable': False, 'renamable': False},
                 {'name': 'l', 'id': 'l',
                  'deletable': False, 'renamable': False},
                 {'name': 'Surface energy (eV/Å^2)', 'id': 'surface_energy',
                  'deletable': False, 'renamable': False}],
        data=[{'h': 1, 'k': 0, 'l': 0, 'surface_energy': 1}],
        editable=True,
        row_deletable=True,
        style_cell={ "textAlign": "center", 'minWidth': '100px'}, 
        fill_width=False
    ),

    # add a button for adding more facets
    html.Button('Add surface', id='editing-rows-button', n_clicks=0),
    dcc.Input(id="MPID", type="text", placeholder="MPID", style={'marginRight':'10px'}),
    dcc.Graph(id='wulff_shape')
])

# app function to allow additional facets upon click
@app.callback(
    Output('hkl_and_surface_energy', 'data'),
    Input('editing-rows-button', 'n_clicks'),
    State('hkl_and_surface_energy', 'data'),
    State('hkl_and_surface_energy', 'columns'))
def add_row(n_clicks, rows, columns):
    if n_clicks > 0:
        rows.append({c['id']: '' for c in columns})
    return rows

# app function to get the lattice parameter and Wulff shape
@app.callback(
    Output('wulff_shape', 'figure'),
    Input('hkl_and_surface_energy', 'data'),
    Input('abc', 'data'),
    Input('angles', 'data'),
    Input("MPID", "value"))
def display_wulff_shape(hkl_and_se, abc, angles, mpid=None):
    
    if mpid:
        surface_data = mpr.get_surface_data(mpid)
        miller_indices = [tuple(surf['miller_index']) for surf in surface_data['surfaces']]
        surface_energies = [surf['surface_energy'] for surf in surface_data['surfaces']]
        latt = mpr.get_structure_by_material_id(mpid, conventional_unit_cell=True).lattice 

    else:
        miller_indices = [(int(row['h']), int(row['k']), int(row['l'])) for row in hkl_and_se]
        surface_energies = [float(row['surface_energy']) for row in hkl_and_se]
        abc = abc[0]
        angles = angles[0]
        latt = Lattice.from_parameters(float(abc['a']), float(abc['b']), float(abc['c']), 
                                       float(angles['alpha']), float(angles['beta']), float(angles['gamma']))
        
    wulff = WulffShape(latt, miller_indices, surface_energies)    
        
    return wulff.get_plotly()

if __name__ == '__main__':
    app.run_server(debug=True)
