# coding: utf-8
# python mainModuleSums.py --subdet all_subdet --event -1 --algo 8towers --particle photons

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

    parser.add_argument("--subdet", default='all_subdet', help="Select subdetector (CEE, CEH, CEH_even, CEH_odd, or all_subdet)") #FIXME
    parser.add_argument("--event", default='5492', help="Select event number or -1 for all events")
    parser.add_argument("--algo", default='8towers', help="Select algorithm (baseline, area_overlap, 8towers)")
    parser.add_argument("--particle", default='photons', help="Select particle type (photons or pions)")

    return parser.parse_args()


def main(subdet, event, particle, algo):
    process = processingMS.Processing() 

    # Method that retrieves events and process the data with the geometry
    data = process.get_data_new(event) 
    print("Dataframe columns",data.columns)

    bin_geojson_filename = '/geojson/towers_bins.geojson'
    hex_geojson_filename = '/geojson/hexagons_byesplit.geojson'
    hdf5_filename = f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/hdf5_files/overlap_data_final_{particle}_{event}.h5'

    initial_kw = {
        'NbinsEta': 20,
        'NbinsPhi': 72,
        'MinPhi': -3.14159,
        'MaxPhi': +3.14159,
        'EtaMin': 1.305,
        'EtaMax': 3.045
    }

    #towers_bins = process.create_bin_df_new(initial_kw)
    process.ModSumToTowers(initial_kw, data , subdet, event, particle, algo, bin_geojson_filename, hex_geojson_filename, hdf5_filename)

    #process.save_bin_geo(towers_bins, output_file)
    #process.save_bin_hex(output_file_hex)

    #process.plot_bins_with_eta_phi_for_phi_90(input_file, output_dir)


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

    

