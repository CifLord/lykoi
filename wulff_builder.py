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
])

@app.callback(
    Output('wulff_shape', 'figure'),
    Output('hkl_and_surface_energy', 'data'),
    Output('abc', 'data'),
    Output('angles', 'data'),
    Output('MPID', 'value'), # returns nothing in order to clear input box
    Output('editing-rows-button', 'n_clicks'), # resets n_clicks to 0 to avoid creating new rows for each input
    Output('calculate_button','disabled'),
    Output('calculate_button', 'n_clicks'), # resets n_clicks to 0 to avoid creating new rows for each input

    Input('hkl_and_surface_energy', 'data'),
    Input('abc', 'data'),
    Input('angles', 'data'),
    Input('wulff_shape', 'figure'),
    Input("MPID", "value"),
    Input('editing-rows-button', 'n_clicks'),
    Input('calculate_button', 'n_clicks'),
    Input('slab_vrun', 'contents'),
    Input('slab_vrun', 'filename'),
    Input('bulk_vrun', 'contents'),
    Input('bulk_vrun', 'filename'),
    Input('hkl_xml', 'data'))
def display_wulff_shape(hkl_and_se, abc, angles, old_wulff_shape, 
                        mpid=None, n_clicks=0, calculate=0, slab_vrun=None, slab_filename=None,
                        bulk_vrun=None, bulk_filename=None, hkl_xml=None):
    
    columns=[{'name': 'h', 'id': 'h', 'deletable': False, 'renamable': False},
             {'name': 'k', 'id': 'k', 'deletable': False, 'renamable': False},
             {'name': 'l', 'id': 'l', 'deletable': False, 'renamable': False},
             {'name': 'Surface energy', 'id': 'surface_energy', 
              'deletable': False, 'renamable': False},
             {'name': 'Area fraction', 'id': 'area_frac', 
              'deletable': False, 'renamable': False}]  
    
    if slab_vrun:
        content_type, content_string = slab_vrun.split(',')
        decoded = base64.b64decode(content_string)
        with NamedTemporaryFile(suffix=slab_filename) as tmp:
            tmp.write(decoded)
            tmp.flush()
            slab_vrun = Vasprun(tmp.name)
        slab_energy = slab_vrun.final_energy
        slab = slab_vrun.final_structure
    if bulk_vrun:
        content_type, content_string = bulk_vrun.split(',')
        decoded = base64.b64decode(content_string)
        with NamedTemporaryFile(suffix=bulk_filename) as tmp:
            tmp.write(decoded)
            tmp.flush()
            bulk_vrun = Vasprun(tmp.name)
        bulk_energy = bulk_vrun.final_energy
        bulk = bulk_vrun.final_structure
    if calculate > 0:
        hkl = (int(hkl_xml[0]['h']), int(hkl_xml[0]['k']), int(hkl_xml[0]['l']))
        slabentry = SlabEntry(slab, slab_energy, hkl)
        calc_surface_energy = slabentry.surface_energy(ComputedStructureEntry(bulk, bulk_energy))
        hkl_and_se.append({'h': hkl[0], 'k': hkl[1], 'l': hkl[-1],
                           'surface_energy': '%.3f' %(calc_surface_energy*EV_PER_ANG2_TO_JOULES_PER_M2)})

    if n_clicks > 0:
        hkl_and_se.append({c['id']: '' for c in columns if c['id'] != 'area_frac'})
    if mpid:
        surface_data = mpr.get_surface_data(mpid)
        miller_indices = [tuple(surf['miller_index']) for surf in surface_data['surfaces']]
        surface_energies = [surf['surface_energy'] for surf in surface_data['surfaces']]
        
        # reset lattice parameter table for this particular mpid
        latt = mpr.get_structure_by_material_id(mpid, conventional_unit_cell=True).lattice 
        abc = [{'a': latt.a, 'b': latt.b, 'c': latt.c}]
        angles = [{'alpha': latt.alpha, 'beta': latt.beta, 'gamma': latt.gamma}]

        # reset the table for this particular mpid
        hkl_and_se=[]
        for i, hkl in enumerate(miller_indices):
            hkl_and_se.append({'h': hkl[0], 'k': hkl[1], 'l': hkl[-1], 'surface_energy': '%.3f' %(surface_energies[i])})
        
    else:
        # only consider rows with appropriate values for h, k, l and surface energy, ignore otherwise 
        miller_indices = [(int(row['h']), int(row['k']), int(row['l'])) for row in hkl_and_se 
                          if all([v != '' and v != None for v in row.values()])]
        surface_energies = [float(row['surface_energy']) for row in hkl_and_se 
                            if all([v != '' and v != None for v in row.values()])]
        latt = Lattice.from_parameters(float(abc[0]['a']), float(abc[0]['b']), float(abc[0]['c']), 
                                       float(angles[0]['alpha']), float(angles[0]['beta']), float(angles[0]['gamma']))
        
    slab_vrun = None if calculate > 0 else slab_vrun
    bulk_vrun = None if calculate > 0 else bulk_vrun
        
    try:
        wulff = WulffShape(latt, miller_indices, surface_energies)
        # add the area fractions
        for i, row in enumerate(hkl_and_se):
            if all([v != '' and v != None for v in row.values()]):
                hkl_and_se[i]['area_frac'] = '%.3f' %(wulff.area_fraction_dict[tuple([int(row['h']), 
                                                                                int(row['k']), 
                                                                                int(row['l'])])])
        return wulff.get_plotly(), hkl_and_se, abc, angles, '', 0, slab_vrun == None or bulk_vrun == None, 0
    except QhullError:
        # If a Wulff shape cannot be enclosed, return the previous Wulff shape
        return old_wulff_shape, hkl_and_se, abc, angles, '', 0, slab_vrun == None or bulk_vrun == None, 0
    
if __name__ == '__main__':
    app.run_server(debug=True)
