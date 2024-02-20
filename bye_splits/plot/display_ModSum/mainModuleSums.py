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

initial_kw = {
    'NbinsEta': 20,
    'NbinsPhi': 24,
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
    #fig = plotMS.plot_single_event(grid_params, data, output_path)
    print("interactive plotting")
    fig= plotMS.plot_single_event_proj(initial_kw, data, df_bin, '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/today_single_event.png' )
    return fig

if __name__ == '__main__':
    args = parsing.parser_display_plotly()
    process = processingMS.Processing() 
    
    #method that retrieves a single event and process the data with the geometry
    data = process.get_data_new("-1") 
    print(data.columns)

    #data_new = process.get_data_new("-1")
    #print("INFO DATAFRAME", data_new.info)
    print("DF", data)
    #method that creates the x/y grid form eta/phi bins 
    df_bin = process.create_bin_df(initial_kw)
    
    #modules to towers algorithms 
    process.ModSumToTowers(initial_kw, data, df_bin, algo='baseline')

    app.run_server(debug=True, port=8051,host="0.0.0.0")

