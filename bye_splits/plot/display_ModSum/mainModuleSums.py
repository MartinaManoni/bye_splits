# coding: utf-8

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

from bye_splits.utils import parsing, common
import plotMS
import processingMS

#20x 24 eta-phi bin for 120 degrees, so for the entire module is 20x72 (360 degrees)
initial_kw = {
    'NbinsEta': 20,
    'NbinsPhi': 72,
    'MinPhi': -3.14159,
    'MaxPhi': +3.14159,
    'EtaMin': 1.5,
    'EtaMax': 3
}

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

    #df_grid = plotMS.create_grid_df(kw)

    app.config.suppress_callback_exceptions = True
    
    print("interactive plotting")
    fig= plotMS.plot_single_event_proj(initial_kw, data, df_bin, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/grid_eta_1.5_3_20x72{algo}_{subdet}_{particle}_{event}.png', dpi=500)
    return fig

if __name__ == '__main__':
    #args = parsing.parser_display_plotly()
    process = processingMS.Processing() 

    subdet = 'all_subdet'  # select between CEE, CEH_even, CEH_odd or all_subdet
    event= '3737' # select between 'event number' or '-1'(all events)
    algo='baseline' # select between baseline, area_overlap, algo3
    particle= 'pions' #select between photons or pions

    #n.b: algo3 for all_subdet not working properly at the moment 

    # Method that retrieves a single event or all events and process the data with the geometry
    data = process.get_data_new(event) 
    print("Dataframe columns",data.columns)
    
    if subdet != 'all_subdet':
        mask = data['group'] == subdet
        data = data[mask]

    data.rename(columns={'wx_center': 'ts_x', 'wy_center': 'ts_y'}, inplace=True) #only because now I am using byesplit geometry
    #print(f"Dataframe {subdet}", data[mask])

    # Method that creates the x/y grid form eta/phi bins - used only for graphical reasons, cause here the bins are rectangular polygons (no arcs)
    df_bin = process.create_bin_df(initial_kw)
    #print("Bin dataframe", df_bin.columns)

    # Modules To Towers algorithms 
    process.ModSumToTowers(initial_kw, data, df_bin, algo, subdet, event, particle)  #poi inserire il subdet qui dentro 

    #plotMS.plot_single_event_proj(initial_kw, data, df_bin, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/grid_eta_1.5_3_20x72{algo}_{subdet}_{particle}_{event}.png')

    #app.run_server(debug=True, port=8051,host="0.0.0.0")

