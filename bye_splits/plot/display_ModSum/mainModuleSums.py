# coding: utf-8
# python mainModuleSums.py --event 5492 --algo baseline --particle photons --subdet CEE

_all_ = [ ]

import os
import sys

parent_dir = os.path.abspath(__file__ + 4 * '/..')
sys.path.insert(0, parent_dir)

from dash import Dash, dcc, html, Input, Output, State, ctx
from dash.exceptions import PreventUpdate
from dash_bootstrap_templates import load_figure_template
import dash_bootstrap_components as dbc
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

import plotMS
import processingMS


def parse_arguments():
    parser = argparse.ArgumentParser(description="Interactive Grid Comparison")

    parser.add_argument("--subdet", default='all_subdet', help="Select subdetector (CEE, CEH, all_subdet)")
    parser.add_argument("--event", default='5492', help="Select event number or -1 for all events")
    parser.add_argument("--algo", default='8towers', help="Select algorithm (baseline, area_overlap, 8towers)")
    parser.add_argument("--particle", default='photons', help="Select particle type (photons or pions)")

    return parser.parse_args()


def main(subdet, event, particle, algo):
    process = processingMS.Processing() 

    # Method that retrieves events and process the data with the geometry
    data = process.get_data_new(event) 

    file_path = f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/photons_manoni/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'

    #process.read_hdf5_structure(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/photons_manoni/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5')
    #process.read_all_block0_values(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/photons_manoni/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5')

    data_gen = process.get_genpart_data(file_path, event)
    ##print("Dataframe columns",data.columns)

    bin_geojson_filename = '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_with_arcs.geojson'
    hex_geojson_filename = '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/hexagons_CMSSW.geojson'
    scint_geojson_filename = '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/scint_modules_geo.geojson'

    hdf5_filename = f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/hdf5_files/overlap_data_{particle}_{event}_shifted.h5'

    #overlap = process.eval_hex_bin_overlap(data, bin_geojson_filename,  hdf5_filename)

    initial_kw = {
        'NbinsEta': 20,
        'NbinsPhi': 72,
        'MinPhi': -3.14159,
        'MaxPhi': +3.14159,
        'EtaMin': 1.305,
        'EtaMax': 3.045
    }

    bins_data, hexagons_data, scint_data = process.read_geojson_files(bin_geojson_filename, hex_geojson_filename, scint_geojson_filename)
    plotMS.plot_full_geom(bins_data, hexagons_data, scint_data, 'plot_layers', plot_type='all')

    #towers_bins = process.create_bin_df_new(initial_kw)
    ##process.ModSumToTowers(initial_kw, data , subdet, event, particle, algo, bin_geojson_filename, hex_geojson_filename, hdf5_filename, data_gen)

    #process.save_bin_geo(towers_bins, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_with_arcs.geojson', f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_only_vertices.geojson')
    #process.save_bin_hex(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/hexagons_CMSSW.geojson')
    #process.save_scint_mod_geo(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/scint_modules_geo.geojson')

    #Dash app 
    #FIXME --> needs to be updated, layer by layer plotting
    '''
    initial_df_grid = processingMS.Processing().create_grid_df(initial_kw)

    app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])

    app.layout = html.Div([
        html.H1("Interactive Grid Comparison"),
        
        dcc.Graph(id='grid-plot'),

        html.Label("Adjust Binning Parameters:"),
        
        dcc.Slider(
            id='phi-slider',
            min=10,
            max=50,
            step=1,
            marks={i: str(i) for i in [10, 20, 30, 40]},
            value=initial_kw['NbinsPhi'],
            tooltip={'placement': 'bottom'},
            updatemode='mouseup'
        ),

        dcc.Slider(
            id='eta-slider',
            min=5,
            max=20,
            step=1,
            marks={i: str(i) for i in [5, 10, 20]},
            value=initial_kw['NbinsEta'],
            tooltip={'placement': 'bottom'},
            updatemode='mouseup'
        )
    ])

    @app.callback(
        Output('grid-plot', 'figure'),
        [Input('phi-slider', 'value'),
        Input('eta-slider', 'value')]
    )
    
    def update_plot(n_phi, n_eta):
        kw = {
            'NbinsEta': n_eta,
            'NbinsPhi': n_phi,
            'MinPhi': -3.14159,
            'MaxPhi': +3.14159,
            'EtaMin': 1.5,
            'EtaMax': 3
        }

        print("interactive plotting")
        fig= plotMS.plot_single_event_proj(initial_kw, data, df_bin, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/grid_eta_1.5_3_20x72{algo}_{subdet}_{particle}_{event}.png', dpi=500)
        return fig'''

        #df_grid = plotMS.create_grid_df(kw)

    #app.config.suppress_callback_exceptions = True
    #app.run_server(debug=True, port=8051, host="0.0.0.0")


if __name__ == '__main__':
    args = parse_arguments()
    main(args.subdet, args.event, args.particle, args.algo)

    

