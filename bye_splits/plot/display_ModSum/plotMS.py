# coding: utf-8

_all_ = []

import os
from pathlib import Path
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)

#from bye_splits.plot.display_plotly import yaml, np, pd, go, dcc

import h5py
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from scipy.spatial.distance import cdist
from scipy.stats import norm

import matplotlib
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PolygonPlt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from shapely.geometry import Polygon, shape
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

import processingMS
#matplotlib.use("TkAgg")

def plot_grid(df_grid, df_bin, kw):
    fig = go.Figure()

    # Sort the DataFrame for continuous lines
    df_grid = df_grid.sort_values(['Eta', 'Phi'])

    # Create separate traces for each line in Eta
    for eta_value in df_grid['Eta'].unique():
        subset = df_grid[df_grid['Eta'] == eta_value]
        fig.add_trace(go.Scatter(
            x=subset['X'],
            y=subset['Y'],
            mode='lines+markers',
            line=dict(color='black'),  # Adjust the color as needed
            marker=dict(size=2),
            text=subset['Phi'],
            hoverinfo='text',
            name=f'Eta={eta_value}'
        ))

    for phi_value in df_grid['Phi'].unique():
        subset = df_grid[df_grid['Phi'] == phi_value]

        fig.add_trace(go.Scatter(
            x=subset['X'],
            y=subset['Y'],
            mode='lines+markers',
            line=dict(color='black'),  # Adjust the color as needed
            marker=dict(size=2),
            text=subset['Eta'],
            hoverinfo='text',
            name=f'Phi={phi_value}'
        ))

    # Plot center_X and center_Y from df_bin
    '''fig.add_trace(go.Scatter(
        x=df_bin['centroid_X'],
        y=df_bin['centroid_Y'],
        mode='markers',
        marker=dict(color='black', size=1),
        name='Center Points'
    )) '''   

    #for vertices_X, vertices_Y in zip(df_bin['bin_vertex_X'], df_bin['bin_vertex_Y']):
    #    vertices_X = np.array(vertices_X)
    #    vertices_Y = np.array(vertices_Y)
    '''vertices_X = np.array(df_bin['bin_vertex_X'][0])
    vertices_Y = np.array(df_bin['bin_vertex_Y'][0])
    fig.add_trace(go.Scatter(
            x=np.append(vertices_X, vertices_X[0]),  # Close the polygon
            y=np.append(vertices_Y, vertices_Y[0]),  # Close the polygon
            mode='markers',
            marker=dict(color='blue', size=4),  # Adjust the color as needed
            name='Polygon Vertices'
        ))'''
    
    for i in range(len(df_bin)):
        vertices_X = df_bin['bin_vertex_X'][i]
        vertices_Y = df_bin['bin_vertex_Y'][i]

        # Add a trace for each polygon
        fig.add_trace(go.Scatter(
            x=np.append(vertices_X, vertices_X[0]),  # Close the polygon
            y=np.append(vertices_Y, vertices_Y[0]),  # Close the polygon
            mode='markers',
            marker=dict(color='black', size=1),  # Adjust the color as needed
            name=f'Polygon {i+1}'  # Add a unique name for each polygon
        ))

    fig.update_layout(
        title='Grid Plot in X and Y',
        xaxis_title='X',
        yaxis_title='Y',
        xaxis=dict(scaleanchor="y", scaleratio=1),
        yaxis=dict(scaleanchor="x", scaleratio=1),
        showlegend=False
    )
    
    return fig

def plot_byesplit_geom(grid_params, hexagon_proj, df_bin, output_path2):
    grid_df = processingMS.Processing().create_grid_df(grid_params)
    test_fig = plot_grid(grid_df, df_bin, grid_params)

    df_bye_sci = pd.read_hdf("/data_CMS/cms/manoni/L1HGCAL/byesplit_geom/geom_bye.hdf5")   
    print("Columns:")
    print(df_bye_sci.columns)   
    
    # Iterate through the rows of aggregated_df
    for index, row in df_bye_sci.iterrows():
        #if row['layer'] <= 28:
        if row['layer'] > 28 and row['layer'] % 2 == 0:    

            hex_x = row['hex_x']
            hex_y = row['hex_y']
            #print("pt", mip_pt)

            hex_x_closed = hex_x + [hex_x[0]]
            hex_y_closed = hex_y + [hex_y[0]]   
            print("arrivato qui 0")
            # Add the hexagon trace to the existing Plotly figure
            test_fig.add_trace(go.Scatter(
                x=hex_x_closed,
                y=hex_y_closed,
                mode='lines',
                line=dict(color='red')
            ))
            print("arrivato qui")
            # Add the hexagon trace to the existing Plotly figure
            test_fig.add_trace(go.Scatter(
                x=df_bye_sci['wx_center'],
                y=df_bye_sci['wy_center'],
                mode='markers',
                marker=dict(color='red', size=4),  # Adjust the color as needed
                name=f'hex_centroid'
            ))
            print("arrivato qui 2")
            '''test_fig.add_trace(go.Scatter(
                x=df_bye_sci['wx_center'],
                y=df_bye_sci['wx_center'],
                mode='text',
                text=[f'U:{wu}, V:{wv}' for wu, wv in zip(df_bye_sci['waferu'], df_bye_sci['waferv'])],
                textposition='middle center',  # Set the position of the text to the middle of each point
                textfont=dict(size=7, color='black'),  # Adjust the font size and color as needed
            ))'''
            print("arrivato qui 3")
    test_fig.write_image(output_path2, width=2400, height=1600)
    return test_fig

def plot_hexagon_with_bins(hexagon_coords, bin_coords, save_path='hexagon_with_bins.png'):
    fig = go.Figure()

    # Plot hexagon
    fig.add_trace(go.Scatter(
        x=hexagon_coords[0],
        y=hexagon_coords[1],
        mode='markers',
        marker=dict(size=8, color='red'),
        name='Hexagon'
    ))

    # Plot bins
    for idx, bin_coord in enumerate(bin_coords):
        fig.add_trace(go.Scatter(
            x=bin_coord[0],
            y=bin_coord[1],
            mode='markers',
            marker=dict(size=1),
            name=f'Bin {idx + 1}'
        ))

    fig.update_layout(
        title='Hexagon with Bins',
        xaxis_title='X',
        yaxis_title='Y',
        xaxis=dict(scaleanchor="y", scaleratio=1),
        yaxis=dict(scaleanchor="x", scaleratio=1),
    )

    # Save as PNG file
    fig.write_image(save_path)
    print(f"Plot saved as {save_path}")

def plot_shifted_modules(silicon_df_proc, shifted_df, df_geom, df_ts, layer_number):
    '''function to plot the centers of module sums for byesplit geometry 
    and CMSSSW, the two geometry have a relative shift'''

    before_shift = silicon_df_proc[(silicon_df_proc['ts_layer'] == layer_number) & (silicon_df_proc['ts_eta']> 0) ]
    after_shift = shifted_df[(shifted_df['ts_layer'] == layer_number)]   
    df_geom_m = df_geom[(df_geom['ts_layer']==layer_number)]
    df_ts = df_ts[(df_ts['ts_layer']==layer_number) & (df_ts['ts_eta']> 0) ]

    plt.figure(figsize=(8, 6))
    plt.scatter(before_shift['wx_center'], before_shift['wy_center'], label='Byesplit geometry (after merge) BEFORE SHIFT', color='red', s=1)
    plt.scatter(after_shift['wx_center'], after_shift['wy_center'], label='Byesplit geometry (after merge) AFTER SHIFT', color='green', s=1)
    #plt.scatter(before_shift['ts_x'], before_shift['ts_y'], label='ntuples data (after merge)', color='black', s=1)
    #plt.scatter(df_ts['ts_x'], df_ts['ts_y'], label='ntuples data (before merge)', color='cyan', s=1)
    #plt.scatter(df_geom_m['wx_center'], df_geom_m['wy_center'], label='CMSSW geom marco (before emerge)', color='green', s=0.5)

    # Plot wafer U and wafer V for wx_center/wy_center (byesplit geometry)
    for idx, row in before_shift.iterrows():
        plt.text(row['wx_center'], row['wy_center']+2, f"U:{row['ts_wu']}, V:{row['ts_wv']}", fontsize=3, color='red')

    for idx, row in after_shift.iterrows():
        plt.text(row['wx_center'], row['wy_center']+2, f"U:{row['ts_wu']}, V:{row['ts_wv']}", fontsize=3, color='green')

    # Plot wafer U and wafer V for ts_x/ts_y  (ntuples data after data/geometry merge)
    #for idx, row in before_shift.iterrows():
        #plt.text(row['ts_x'], row['ts_y'], f"U={row['ts_wu']}, V={row['ts_wv']}", fontsize=3)

    # Plot wafer U and wafer V for wx_center/wy_center (CMSSW geometry)
    #for idx, row in df_geom_m.iterrows():
        #plt.text(row['wx_center'], row['wy_center']+2, f"U={row['ts_wu']}, V={row['ts_wv']}", fontsize=3, color='green')

    # Plot wafer U and wafer V for ts_x/ts_y (ntuples data before data/geometry merge)
    #for idx, row in df_ts.iterrows():
        #plt.text(row['ts_x'], row['ts_y'], f"U:{row['ts_wu']}, V:{row['ts_wv']}", fontsize=3, color='cyan')

    # Plot the hexagons positions before the shift
    for hex_x, hex_y in zip(before_shift['hex_x'], before_shift['hex_y']):
        # Closing the loop by connecting back to the first vertex
        hex_x.append(hex_x[0])
        hex_y.append(hex_y[0])
        plt.plot(hex_x, hex_y, color='red', alpha=0.5)

    # Plot the hexagons positions after the shift
    for hex_x, hex_y in zip(after_shift['hex_x'], after_shift['hex_y']):
        # Closing the loop by connecting back to the first vertex
        hex_x.append(hex_x[0])
        hex_y.append(hex_y[0])
        plt.plot(hex_x, hex_y, color='black', alpha=0.5) 

    plt.title(f'Hexagons positions before and after shift (Layer {layer_number} eta >0)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/Layer{layer_number}_NEW_MART_2.png", dpi = 400)

def plot_layers(df_ts, ds_new):
    #layers = set(df_ts[df_ts['ts_layer'] > 37]['ts_layer'].unique()).intersection(set(ds_new[ds_new['ts_layer'] > 37]['ts_layer'].unique()))
    layers = set(df_ts['ts_layer'].unique()).intersection(set(ds_new['ts_layer'].unique()))
    print(layers)
    for layer in layers:
        fig, ax = plt.subplots()

        # Plot scatter plot for ds_new for the current layer
        ds_new_layer = ds_new[ds_new['ts_layer'] == layer]
        ax.scatter(ds_new_layer['wx_center'], ds_new_layer['wy_center'], color='blue', s=1, label='byesplit geom')
        for index, row in ds_new_layer.iterrows():
            ax.annotate(f"U:{row['ts_wu']}, V:{row['ts_wv']}", (row['wx_center'], row['wy_center']), fontsize=3, color='blue')

        # Plot scatter plot for df_ts for the current layer
        df_ts_layer = df_ts[df_ts['ts_layer'] == layer]
        ax.scatter(df_ts_layer['ts_x'], df_ts_layer['ts_y'], color='red', s=1, label='nuples data')
        for index, row in df_ts_layer.iterrows():
            ax.annotate(f"U:{row['ts_wu']:.0f}, V:{row['ts_wv']:.0f}, E:{row['ts_pt']:.2f}", (row['ts_x'], row['ts_y']+4), fontsize=3, color='red')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'Scatter plot of layer {layer}')
        ax.legend()
        plt.savefig(f"Scatter_plot_of_layer_{layer}.png", dpi=400)
        plt.close()

def plot_layers_sci(df_tc, sci_geom):
    #layers = set(df_ts[df_ts['ts_layer'] > 37]['ts_layer'].unique()).intersection(set(ds_new[ds_new['ts_layer'] > 37]['ts_layer'].unique()))
    #layers = set(df_tc['tc_layer'].unique()).intersection(set(sci_geom['tc_layer'].unique()))
    #print(layers)
    for layer in range(37,51):
        fig, ax = plt.subplots()

        # Plot scatter plot for ds_new for the current layer
        geom_layer = sci_geom[sci_geom['tc_layer'] == layer]
        ax.scatter(geom_layer['x'], geom_layer['y'], color='blue', s=1, label='byesplit geom')
        for index, row in geom_layer.iterrows():
            ax.annotate(f"U:{row['tc_cu']}, V:{row['tc_cv']}", (row['x'], row['y']), fontsize=1, color='blue')

        # Plot scatter plot for df_tc for the current layer
        tc_layer = df_tc[(df_tc['tc_layer'] == layer)] #& (df_tc['tc_cv'] >10 )
        ax.scatter(tc_layer['tc_x'], tc_layer['tc_y'], color='red', s=1, label='nuples data')
        for index, row in tc_layer.iterrows():
            ax.annotate(f"U:{row['tc_cu']:.0f}, V:{row['tc_cv']:.0f}, E:{row['tc_pt']:.2f}", (row['tc_x'], row['tc_y']+4), fontsize=3, color='red')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'TC - Scatter plot of layer {layer}')
        ax.legend()
        plt.savefig(f"sci_plot_of_layer_{layer}.png", dpi=600)
        plt.close()


def plot_baseline(df_baseline_proj, algo, event, particle):
    fig, ax = plt.subplots(figsize=(10, 8))

    # Dictionary to store total pt values for each bin
    total_pt = {}

    for index, row in df_baseline_proj.iterrows():
        eta_vertices = row['eta_vertices']
        phi_vertices = row['phi_vertices']
        pt = row['pt']

        # Store total pt for the bin
        bin_key = (eta_vertices, phi_vertices)
        total_pt[bin_key] = pt

    # Normalize total pt values for color mapping
    max_pt = max(total_pt.values())
    norm = Normalize(vmin=0, vmax=max_pt)

    # Plot bins with total pt values
    for bin_key, pt in total_pt.items():
        eta_vertices, phi_vertices = bin_key
        poly = PolygonPlt(np.column_stack((eta_vertices, phi_vertices)), closed=True, edgecolor='black')
        ax.add_patch(poly)

        # Set color for the bin based on total pt value
        color = ScalarMappable(norm=norm, cmap='viridis').to_rgba(pt)
        poly.set_facecolor(color)

        # Add text annotation with total pt value
        text_x = np.mean(eta_vertices)
        text_y = np.mean(phi_vertices)
        ax.text(text_x, text_y, f'{pt:.1f}', color='black', ha='center', va='center', fontsize=5)

    # Add color bar
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    sm = ScalarMappable(norm=norm, cmap='viridis')
    sm.set_array(list(total_pt.values()))
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Total pt')

    # Set labels and title
    ax.set_xlabel('Eta')
    ax.set_ylabel('Phi')
    ax.set_title(f'{algo}_{particle}_{event}')

    ax.autoscale()
    plt.savefig(f'{algo}_{particle}_{event}.png', dpi=500)  # Save the plot as an image
    plt.show()

def plot_towers_eta_phi_grid(df_baseline_proj, data_gen, algo, event, particle, subdet, results_df):
    #this works only for a single event : need to put options to gen data to take only the n events that coprrespond to one from data
    print("plotting eta_phi_towers")
    fig, ax = plt.subplots(figsize=(10, 8))

    #print ("DATA GEN", data_gen)
    #print ("DATA ", df_baseline_proj)
    # Plotting the grid of bins
    initial_kw = {
        'NbinsEta': 20,
        'NbinsPhi': 72,
        'MinPhi': -3.14159,
        'MaxPhi': +3.14159,
        'EtaMin': 1.305,
        'EtaMax': 3.045
    }
    eta_bins = np.linspace(initial_kw['EtaMin'], initial_kw['EtaMax'], initial_kw['NbinsEta'] + 1)
    phi_bins = np.linspace(initial_kw['MinPhi'], initial_kw['MaxPhi'], initial_kw['NbinsPhi'] + 1)

    for i in range(len(eta_bins) - 1):
        for j in range(len(phi_bins) - 1):
            eta_vertices = [eta_bins[i], eta_bins[i + 1], eta_bins[i + 1], eta_bins[i]]
            phi_vertices = [phi_bins[j], phi_bins[j], phi_bins[j + 1], phi_bins[j + 1]]
            poly = PolygonPlt(np.column_stack((eta_vertices, phi_vertices)), closed=True, edgecolor='black')#, closed=True
            color = ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='viridis').to_rgba(0)
            poly.set_facecolor(color)
            #ax.add_patch(poly)

    # Normalize total pt values for color mapping
    max_pt = df_baseline_proj['pt'].max()
    norm = Normalize(vmin=0, vmax=max_pt)

    # Plotting bins with colors and annotations
    for _, row in df_baseline_proj.iterrows():

        eta_vertices = row['eta_vertices']
        phi_vertices = row['phi_vertices']
        pt = row['pt']

        bin_key = (eta_vertices, phi_vertices)
        color = ScalarMappable(norm=norm, cmap='viridis').to_rgba(pt)

        poly = PolygonPlt(np.column_stack((eta_vertices, phi_vertices)), closed=True, edgecolor='black')#closed=True
        ax.add_patch(poly)
        poly.set_facecolor(color)

        text_x = np.mean(eta_vertices)
        text_y = np.mean(phi_vertices)
        ax.text(text_x, text_y, f'{pt:.1f}', color='black', ha='center', va='center', fontsize=5)
        #ax.text(data_gen['gen_eta'], data_gen['gen_phi'], f'.', color='red', ha='center', va='center', fontsize=5)
        #plt.scatter(data_gen['gen_eta'], data_gen['gen_phi'], marker='o', color='red')
    # Plot gen data points (handle multiple points if data_gen has more than one row)
    if len(data_gen) > 1:
        for idx, row in data_gen.iterrows():
            ax.text(row['gen_eta'], row['gen_phi'], '.', color='red', ha='center', va='center', fontsize=8)
    else:
        ax.text(data_gen['gen_eta'].values[0], data_gen['gen_phi'].values[0], '.', color='red', ha='center', va='center', fontsize=8)

    if particle == "pions":
        # Plot jets as circles based on results_df
        for _, jet_row in results_df.iterrows():
            jet_eta = jet_row['reco_eta']
            jet_phi = jet_row['reco_phi']
            print("jet_eta", jet_eta)
            print("jet_phi", jet_phi)
            #jet_radius = 0.4  # This is the jet radius (Anti-kt algorithm radius)

            # Create a circle for each jet
            circle = Circle((jet_eta, jet_phi), 0.4 , color='red', fill=False, linewidth=2, label=r'Jet (Anti-kt), $\Delta R = 0.4$')
            ax.add_patch(circle)
            
            # Create a circle for each jet
            circle2 = Circle((jet_eta, jet_phi), 0.8 , color='blue', fill=False, linewidth=2, label=r'Jet (Anti-kt), $\Delta R = 0.8$' )
            ax.add_patch(circle2)

            # Create a circle for each jet - FOR MATCHING
            circle2 = Circle((jet_eta, jet_phi), 0.1 , color='green', fill=False, linewidth=2, label=r'Jet (Anti-kt), $\Delta R = 0.1$' )
            ax.add_patch(circle2)

    # Avoid duplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=10)

    # Add color bar
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    sm = ScalarMappable(norm=norm, cmap='viridis')
    sm.set_array(df_baseline_proj['pt'])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Total pt')

    # Set labels and title
    #ax.set_xlabel('Eta')
    #ax.set_ylabel('Phi')
    ax.set_xlabel(r'$\eta$', fontsize=14)
    ax.set_ylabel(r'$\phi$', fontsize=14)
    ax.set_title(f'{algo}_{particle}_{event}')

    ax.autoscale()
    #plt.legend()
    plt.savefig(f'{algo}_{particle}_{event}_{subdet}_eta_phi_towers.png', dpi=500)  # Save the plot as an image
    plt.show()



def plot_hex_bin_overlap_save(hdf5_filename, output_folder):
        with h5py.File(hdf5_filename, 'r') as hf:
            for layer_name in hf.keys():
                layer_group = hf[layer_name]
                plt.figure(figsize=(8, 6))
                plt.title(f'Layer: {layer_name}')
                
                for hex_key in layer_group.keys():
                    hex_group = layer_group[hex_key]
                    hex_x = hex_group['hex_x'][:]
                    hex_y = hex_group['hex_y'][:]
                    
                    # Repeat the first vertex to close the polygon
                    hex_x_closed = list(hex_x) + [hex_x[0]]
                    hex_y_closed = list(hex_y) + [hex_y[0]]
                    
                    plt.plot(hex_x_closed, hex_y_closed, color='blue')

                    # Retrieve pt value from the hexagon's group
                    pt = hex_group['ts_pt'][()]
                    
                    # Retrieve centroid coordinates from the hexagon's group
                    hex_center_x = hex_group['hex_x_centroid'][()]
                    hex_center_y = hex_group['hex_y_centroid'][()]
                    
                    plt.text(hex_center_x, hex_center_y, f'pt: {pt:.2f}', color='red', ha='center', va='center')
                    
                    for bin_key in hex_group.keys():
                        if bin_key.startswith('bin_'):
                            bin_group = hex_group[bin_key]
                            x_vertices = bin_group['x_vertices'][:]
                            y_vertices = bin_group['y_vertices'][:]

                            # Calculate area overlap percentage
                            area_overlap_percentage = bin_group['percentage_overlap'][()]
                            
                            plt.fill(x_vertices, y_vertices, color='orange', alpha=0.5)

                            # Plot area overlap percentage near the bin
                            bin_center_x = np.mean(x_vertices)
                            bin_center_y = np.mean(y_vertices)
                            #plt.text(bin_center_x, bin_center_y, f'{area_overlap_percentage:.2f}', color='black', ha='center', va='center')
                
                plt.xlabel('X')
                plt.ylabel('Y')
                plt.gca().set_aspect('equal', adjustable='box')
                plt.grid(True)
                
                output_filename = os.path.join(output_folder, f'{layer_name}_overlap_2.png')
                plt.savefig(output_filename)
                plt.close()    

def plot_scint_tiles(df):
    """
    This function plots from the original ds_geom['sci'] geometry from Bruno all the TCs of the scintillator
    """
    for layer in df['tc_layer'].unique():
        print("layer num", layer)
        layer_df = df[df['tc_layer'] == layer]

        # Extracting x, y, diamond_x, and diamond_y for the current layer
        x = layer_df['x']
        y = layer_df['y']
        eta = layer_df['eta']
        phi = layer_df['phi']
        ts_ieta = layer_df['ts_ieta']
        ts_iphi = layer_df['ts_iphi']
        tc_cu = layer_df['tc_cu']
        tc_cv = layer_df['tc_cv']
        diamond_x = layer_df['diamond_x']
        diamond_y = layer_df['diamond_y']

        # Plotting x and y variables
        plt.plot(x, y, 'o', label='x, y',markersize=0.5)

        for i in range(len(x)):
            #plt.text(x.iloc[i], y.iloc[i], f'{tc_cu.iloc[i]},{tc_cv.iloc[i]}', ha='center', va='center', fontsize=1, color='black')
            #plt.text(x.iloc[i], y.iloc[i], f'{eta.iloc[i]:.2f},{phi.iloc[i]:.2f}', ha='center', va='center', fontsize=1, color='black')
            plt.text(x.iloc[i], y.iloc[i], f'{ts_ieta.iloc[i]}, {ts_iphi.iloc[i]}', ha='center', va='center', fontsize=1, color='black')

        # Plotting diamond_x and diamond_y as polygons
        for i in range(len(diamond_x)):
            plt.fill(diamond_x.iloc[i], diamond_y.iloc[i], alpha=0.5)

        # Adding labels and legend
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'TC Layer {layer}')
        plt.legend()
        plt.grid(True)
        plt.axis('equal')

        # If filename is provided, save the plot as PNG
        output_dir = '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/plot_scint/'
        output_file = os.path.join(output_dir, f'Layer_{layer}_scintillator_ieta_iphi.png')
        plt.savefig(output_file, dpi=700)
        plt.close()
        print(f"Plot for layer {layer} saved as {output_file}")

def plot_full_geom(bins_data, hexagons_data, scint_data, output_dir, plot_type='all'):
    print("Plotting full geometry")
    layer_data = defaultdict(lambda: {'bins': [], 'hexagons': [], 'scint_mod': []})

    for geojson_data, key in zip([bins_data, hexagons_data, scint_data], ['bins', 'hexagons', 'scint_mod']):
        for feature in geojson_data['features']:
            layer = feature['properties']['Layer']
            geometry = shape(feature['geometry'])
            properties = feature['properties']
            layer_data[layer][key].append({'geometry': geometry, **properties})

    for layer, data in layer_data.items():
        fig, ax = plt.subplots(figsize=(8, 6))

        if plot_type == 'all' or plot_type == 'bins':
            for i, bin_info in enumerate(data['bins']):
                bin_geometry = bin_info['geometry']
                eta_vertices = bin_info['Eta_vertices']
                x, y = bin_geometry.exterior.xy
                ax.plot(x, y, color='black', linewidth=0.2)
                if i < 20:
                    bottom_right_eta = round(eta_vertices[3], 2)
                    bottom_right_x, bottom_right_y = bin_geometry.exterior.coords[3]
                    ax.text(bottom_right_x, bottom_right_y, f'{bottom_right_eta}', ha='center', va='top', fontsize=2, color='black',  rotation=+90)

        if plot_type == 'all' or plot_type == 'hexagons':
            for info in data['hexagons']:
                geometry = info['geometry']
                x, y = geometry.exterior.xy
                ax.plot(x, y, color='red', linewidth=0.7)
                centroid = geometry.centroid
                ax.text(centroid.x, centroid.y, f"{info['wu']},{info['wv']}", ha='center', va='center', fontsize=4, color='red')

        if plot_type == 'all' or plot_type == 'scint_mod':
            for info in data['scint_mod']:
                geometry = info['geometry']
                x, y = geometry.exterior.xy
                ax.plot(x, y, color='blue', linewidth=0.7)
                centroid = geometry.centroid
                ax.text(centroid.x, centroid.y, f"{info['ieta']},{info['iphi']}", ha='center', va='center', fontsize=4, color='blue')

        ax.set_title(f'Layer {layer} - Full geometry')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.grid(True)
        ax.axis('equal')

        output_file = os.path.join(output_dir, f'layer_{layer}_full_geometry_V11_{plot_type}.png')
        fig.savefig(output_file, dpi=700)
        plt.close(fig)

def plot_window_with_subwindows(window_bins, eta_min, eta_max, phi_min, phi_max, eta_start, eta_end, phi_start, phi_end,particle_eta,particle_phi):
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot the 12x12 window
    for idx, row in window_bins.iterrows():
        polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='black', facecolor='none')
        ax.add_patch(polygon)

    # Plot the particle position with a red cross
    plt.scatter(particle_eta, particle_phi, color='red', marker='x')


    # Display the plot
    #plt.xlim(eta_min, eta_max)
    #plt.ylim(phi_min, phi_max)
    plt.xlabel('Eta')
    plt.ylabel('Phi')
    plt.title(f'3x3 Subwindow from ({eta_start}, {phi_start}) to ({eta_end}, {phi_end})')
    plt.savefig(f'Subwindow_from_{eta_start}_{phi_start}_{eta_end}_{phi_end}.png')
    plt.show()
    plt.close()


def plot_window_with_wraparound_only_window(window_bins_part1, window_bins_part2, eta_min, eta_max, phi_min, phi_max, particle_eta, particle_phi):
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot the Part 1 bins (from phi_min to +pi)
    for idx, row in window_bins_part1.iterrows():
        polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='blue', facecolor='none', alpha=0.7, label='Part 1' if idx == 0 else "")
        ax.add_patch(polygon)

    # Plot the Part 2 bins (from -pi to phi_max)
    for idx, row in window_bins_part2.iterrows():
        polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='red', facecolor='none', alpha=0.7, label='Part 2' if idx == 0 else "")
        ax.add_patch(polygon)

    # Plot the overall eta and phi limits
    ax.axhline(y=np.pi, color='grey', linestyle='--', label='Phi = +pi')
    ax.axhline(y=-np.pi, color='grey', linestyle='--', label='Phi = -pi')

    # Highlight the eta and phi boundaries
    ax.plot([eta_min, eta_max, eta_max, eta_min, eta_min], [phi_min, phi_min, phi_max, phi_max, phi_min], 'k--', label='Window Bounds')

    # Plot the generated particle's eta and phi position
    ax.scatter(particle_eta, particle_phi, color='black', marker='x', s=100, label='Generated Particle')

    # Set labels and title
    ax.set_xlabel('Eta Center')
    ax.set_ylabel('Phi Center')
    ax.set_title('Window Bins with Wrap-Around in Phi')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)

    # Set axis limits for better visualization
    ax.set_xlim(eta_min - 0.1, eta_max + 0.1)
    ax.set_ylim(-np.pi - 0.1, np.pi + 0.1)

    plt.show()

def plot_window_with_wraparound(window_bins, window_bins_part1, window_bins_part2,  eta_min, eta_max, phi_min, phi_max, particle_eta_1, particle_phi_1, eta_start, eta_end, phi_start, phi_end):
    """
    Plots the main window (12x12) and a moving subwindow (3x3) to visualize how the subwindows are iterated over.

    Parameters:
    - window_bins_part1: DataFrame for the first part of the window, handling wrap-around cases
    - window_bins_part2: DataFrame for the second part of the window, handling wrap-around cases
    - eta_min, eta_max: Bounds of the main window in eta
    - phi_min, phi_max: Bounds of the main window in phi
    - particle_eta, particle_phi: Coordinates of the particle to plot as a red cross
    - eta_start, eta_end: Bounds of the subwindow in eta
    - phi_start, phi_end: Bounds of the subwindow in phi
    """

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 10))

    if window_bins is not None:
        # Plot the main window bins (part 1)
        for idx, row in window_bins.iterrows():
            polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='black', facecolor='none')
            ax.add_patch(polygon)
    else:
        if window_bins_part1 is not None:
            # Plot the main window bins (part 1)
            for idx, row in window_bins_part1.iterrows():
                polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='black', facecolor='none')
                ax.add_patch(polygon)
                
        if window_bins_part2 is not None:
            # Plot the main window bins (part 2)
            for idx, row in window_bins_part2.iterrows():
                polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='black', facecolor='none')
                ax.add_patch(polygon)

    # Highlight the subwindow in red
    if phi_start > phi_end:
        # Handle wrap-around case for the subwindow
        # Part 1: From phi_start to +pi
        subwindow_bins_part1 = window_bins_part1[
            (window_bins_part1['eta_center'] >= eta_start) & (window_bins_part1['eta_center'] <= eta_end) &
            (window_bins_part1['phi_center'] >= phi_start) & (window_bins_part1['phi_center'] <= np.pi)
        ]
        subwindow_bins_part2 = window_bins_part2[
            (window_bins_part2['eta_center'] >= eta_start) & (window_bins_part2['eta_center'] <= eta_end) &
            (window_bins_part2['phi_center'] >= -np.pi) & (window_bins_part2['phi_center'] <= phi_end)
        ]
        for idx, row in subwindow_bins_part1.iterrows():
            polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='red', facecolor='none')
            ax.add_patch(polygon)
        for idx, row in subwindow_bins_part2.iterrows():
            polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='red', facecolor='none')
            ax.add_patch(polygon)
    else:
        # Normal case for the subwindow
        subwindow_bins = window_bins[
            (window_bins['eta_center'] >= eta_start) & (window_bins['eta_center'] <= eta_end) &
            (window_bins['phi_center'] >= phi_start) & (window_bins['phi_center'] <= phi_end)
        ]
        for idx, row in subwindow_bins.iterrows():
            polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='red', facecolor='none')
            ax.add_patch(polygon)

    # Plot the particle position with a red cross
    ax.scatter(particle_eta_1, particle_phi_1, color='red', marker='x', label='Particle Position')

    # Set axis labels and title
    ax.set_xlabel('Eta')
    ax.set_ylabel('Phi')
    ax.set_title('Window and Subwindow Visualization')

    # Set axis limits
    # Set axis limits for better visualization
    ax.set_xlim(eta_min - 0.1, eta_max + 0.1)
    ax.set_ylim(-np.pi - 0.1, np.pi + 0.1)

    # Add legend
    ax.legend()

    # Show the plot
    plt.close()
    #plt.savefig("prova_window.png")


def plot_eta_phi_resolution(results_df, algo, event, particle, subdet):
    """
    Plot histograms of eta and phi differences for each event with Gaussian fits.
    
    Parameters:
    - results_df: DataFrame containing 'eta_diff' and 'phi_diff' columns for each event.
    """
    eta_diffs = results_df['eta_diff']
    phi_diffs = results_df['phi_diff']

    # Plot histogram for eta differences
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    n, bins, patches = plt.hist(eta_diffs, bins=10, color='blue', alpha=0.4, density=True)
    mu_eta, std_eta = norm.fit(eta_diffs)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu_eta, std_eta)
    plt.plot(x, p, 'k', linewidth=2)
    plt.xlabel('Eta Difference')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Eta Differences\n mu={mu_eta:.5f}, sigma={std_eta:.5f}')
    #plt.savefig('eta_diff_histogram.png')

    # Plot histogram for phi differences
    plt.subplot(1, 2, 2)
    n, bins, patches = plt.hist(phi_diffs, bins=10, color='green', alpha=0.4, density=True)
    mu_phi, std_phi = norm.fit(phi_diffs)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu_phi, std_phi)
    plt.plot(x, p, 'k', linewidth=2)
    plt.xlabel('Phi Difference')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Phi Differences\n mu={mu_phi:.5f}, sigma={std_phi:.5f}')
    plt.savefig(f'{algo}_{particle}_{event}_{subdet}_eta_phi_resolution_8x8.png')

    plt.tight_layout()
    plt.show()

def plot_energy_ratio_histogram():
        # Read the energy ratios from the file
        pt_ratios = []
        with open('pt_ratio_overlap.txt', 'r') as file:
            for line in file:
                try:
                    ratio = float(line.strip())
                    pt_ratios.append(ratio)
                except ValueError:
                    # Handle the case where the line is not a valid float
                    continue

        # Plot the histogram
        plt.figure(figsize=(10, 6))
        plt.hist(pt_ratios, bins=50, edgecolor='black', alpha=0.7)
        plt.title('Histogram of pt Ratios')
        plt.xlabel('pt Ratio')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.show()

#checks V16 integration
def plot_layer_hexagons(df, layer, x1,y1):
    group = df[df['ts_layer'] == layer]
    if group.empty:
        print(f"No data for layer {layer}")
        return

    fig, ax = plt.subplots()
    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(group['ts_pt'].min(), group['ts_pt'].max())

    for _, row in group.iterrows():
        hex_x, hex_y = row['hex_x'], row['hex_y']
        polygon = Polygon(zip(hex_x, hex_y))
        color = cmap(norm(row['ts_pt']))
        x, y = polygon.exterior.xy
        ax.fill(x, y, color=color)

    ax.plot(x1, y1,'rx')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label='ts_pt')
    ax.set_title(f'Layer {layer}')
    print(x1, y1)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()


def plot_towers_xy_grid(df_baseline_proj, data_gen, algo, event, particle, subdet):
    print("plotting xy_towers")
    fig, ax = plt.subplots(figsize=(10, 8))

    #print("DATA GEN", data_gen)
    # Plotting the grid of bins
    initial_kw = {
        'NbinsEta': 20,
        'NbinsPhi': 72,
        'MinPhi': -3.14159,
        'MaxPhi': +3.14159,
        'EtaMin': 1.305,
        'EtaMax': 3.045
    }
    eta_bins = np.linspace(initial_kw['EtaMin'], initial_kw['EtaMax'], initial_kw['NbinsEta'] + 1)
    phi_bins = np.linspace(initial_kw['MinPhi'], initial_kw['MaxPhi'], initial_kw['NbinsPhi'] + 1)

    # Normalize total pt values for color mapping
    max_pt = df_baseline_proj['pt'].max()
    norm = Normalize(vmin=0, vmax=max_pt)

    # Plotting bins with colors and annotations
    for _, row in df_baseline_proj.iterrows():
        eta_vertices = row['eta_vertices']
        phi_vertices = row['phi_vertices']
        pt = row['pt']

        x_vertices, y_vertices = [], []
        for eta, phi in zip(eta_vertices, phi_vertices):
            x, y = sph2cart(eta, phi)
            x_vertices.append(x)
            y_vertices.append(y)

       #color = ScalarMappable(norm=norm, cmap='viridis').to_rgba(pt)
        if pt > 0:
            color = ScalarMappable(norm=norm, cmap='viridis').to_rgba(pt)
        else:
            color = 'none'
        poly = PolygonPlt(np.column_stack((x_vertices, y_vertices)), closed=True, edgecolor='black')#edgecolor='black'
        poly.set_facecolor(color)
        ax.add_patch(poly)

        text_x = np.mean(x_vertices)
        text_y = np.mean(y_vertices)
        #ax.text(text_x, text_y, f'{pt:.1f}', color='black', ha='center', va='center', fontsize=5)

    # Plot data_gen point
    for _, row in data_gen.iterrows():
        gen_x, gen_y = sph2cart(row['gen_eta'], row['gen_phi'])
        ax.text(gen_x, gen_y, 'x', color='red', ha='center', va='center', fontsize=9, fontweight='bold')
    #gen_x, gen_y = sph2cart(data_gen['gen_eta'], data_gen['gen_phi'])
    #ax.text(gen_x, gen_y, f'x', color='red', ha='center', va='center', fontsize=9, fontweight='bold')

    # Add color bar
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    sm = ScalarMappable(norm=norm, cmap='viridis')
    sm.set_array(df_baseline_proj['pt'])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Total pT [GeV]')

    # Set labels and title
    ax.set_xlabel('X poistion [cm]')
    ax.set_ylabel('Y position [cm]')
    ax.set_title(f'{algo}_{particle}_{event}')

    ax.set_aspect('equal')
    ax.autoscale()
    plt.savefig(f'{algo}_{particle}_{event}_{subdet}_xy_towers.png', dpi=700)  # Save the plot as an image
    plt.show()


def plot_hex_bins_0(df_hexagon_info):
        # Get the unique events from the index
        unique_events = df_hexagon_info.index.get_level_values('event').unique()

        # Select the first event (or specify another event if needed)
        first_event_id = unique_events[0]
        print(f"Plotting for event: {first_event_id}")

        # Filter the DataFrame for the selected event
        event_data = df_hexagon_info.xs(first_event_id, level='event')

        # Iterate over layers in the selected event
        for layer, layer_data in event_data.groupby('layer'):
            plt.figure(figsize=(10, 10))
            plt.title(f"Event {first_event_id}, Layer {layer}")
            plt.xlabel("X Coordinates")
            plt.ylabel("Y Coordinates")

            # Iterate over hexagons in the layer and plot them
            for _, hexagon in layer_data.iterrows():
                hex_x = hexagon['hex_x']
                #print("hex_x", hex_x)
                hex_y = hexagon['hex_y']
                #print("hex_y", hex_y)
                polygon = Polygon(zip(hex_x, hex_y))
                x, y = polygon.exterior.xy  # Get the exterior coordinates of the polygon

                # Plot the hexagon
                plt.fill(x, y, alpha=0.5, edgecolor='black')

                # Plot overlapping bins, if any
                for bin_info in hexagon['bins_overlapping']:
                    bin_x = bin_info['x_vertices']
                    bin_y = bin_info['y_vertices']
                    plt.fill(bin_x, bin_y, alpha=0.3, color='red', edgecolor='black')

                    # Calculate the centroid of the bin for placing the text annotation
                    bin_polygon = Polygon(zip(bin_x, bin_y))
                    centroid_x, centroid_y = bin_polygon.centroid.x, bin_polygon.centroid.y

                    # Annotate the percentage overlap
                    percentage_overlap = bin_info['percentage_overlap']
                    plt.text(centroid_x, centroid_y, f"{percentage_overlap:.2f}%",
                            fontsize=9, ha='center', va='center', color='blue')

            plt.savefig(f"event_{first_event_id}_layer_{layer}_subdet1.png", format='png')
            plt.close()  # Close the plot to free up memory and avoid displaying it

def plot_hex_bins(df_hexagon_info):
    print("df_hexagon_info", df_hexagon_info.columns)
    # Get the unique events from the index
    unique_events = df_hexagon_info.index.get_level_values('event').unique()

    # Select the first event (or specify another event if needed)
    first_event_id = unique_events[0]
    print(f"Plotting for event: {first_event_id}")

    # Filter the DataFrame for the selected event
    event_data = df_hexagon_info.xs(first_event_id, level='event')

    # Iterate over layers in the selected event
    for layer, layer_data in event_data.groupby('layer'):
        plt.figure(figsize=(10, 10))
        plt.title(f"Event {first_event_id}, Layer {layer}")
        plt.xlabel("X Coordinates")
        plt.ylabel("Y Coordinates")

        # Iterate over hexagons in the layer and plot them
        for _, hexagon in layer_data.iterrows():
            hex_x = hexagon['hex_x']
            hex_y = hexagon['hex_y']
            polygon = Polygon(zip(hex_x, hex_y))
            x, y = polygon.exterior.xy  # Get the exterior coordinates of the polygon

            # Plot the hexagon
            plt.fill(x, y, alpha=0.5, edgecolor='black')

            # Calculate the centroid of the hexagon
            centroid_x, centroid_y = polygon.centroid.x, polygon.centroid.y

            # Convert centroid to (eta, phi) spherical coordinates
            eta, phi = cart2sph(centroid_x, centroid_y, 350.)

            # Annotate the centroid with eta and phi
            plt.text(centroid_x, centroid_y, f"η={eta:.2f}, φ={phi:.2f}",
                     fontsize=9, ha='center', va='center', color='black')

            # Plot overlapping bins, if any
            for bin_info in hexagon['bins_overlapping']:
                bin_x = bin_info['x_vertices']
                bin_y = bin_info['y_vertices']
                plt.fill(bin_x, bin_y, alpha=0.3, color='red', edgecolor='black')

                # Calculate the centroid of the bin for placing the text annotation
                bin_polygon = Polygon(zip(bin_x, bin_y))
                bin_centroid_x, bin_centroid_y = bin_polygon.centroid.x, bin_polygon.centroid.y
                eta_bin, phi_bin= cart2sph(bin_centroid_x, bin_centroid_y, 350.)

                # Annotate the percentage overlap
                percentage_overlap = bin_info['percentage_overlap']
                plt.text(bin_centroid_x, bin_centroid_y, f"{percentage_overlap:.2f}%",
                         fontsize=9, ha='center', va='center', color='blue')
                plt.text(bin_centroid_x, bin_centroid_y, f"η={eta_bin:.2f}, φ={phi_bin:.2f}",
                         fontsize=9, ha='center', va='center', color='black')

        plt.savefig(f"event_{first_event_id}_layer_{layer}_subdet1.png", format='png')
        plt.close()  # Close the plot to free up memory and avoid displaying it


def sph2cart(eta, phi, z=322.):
        ''' Useful conversion: Spherical coordinates to cartesian coordinates (x, y)  '''
        theta = 2*np.arctan(np.exp(-eta))
        r = z / np.cos(theta)
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        return x, y

def cart2sph(x, y, z=322.):
        ''' Useful conversion: Cartesian coordinates to spherical coordinates (eta, phi) '''
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z / r)
        eta = -np.log(np.tan(theta / 2))
        phi = np.arctan2(y, x)
        return eta, phi