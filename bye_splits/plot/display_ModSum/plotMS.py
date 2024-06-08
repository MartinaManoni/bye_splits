# coding: utf-8

_all_ = []

import os
from pathlib import Path
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)

from bye_splits.plot.display_plotly import yaml, np, pd, go, dcc

import h5py
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from scipy.spatial.distance import cdist
from scipy.stats import norm

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PolygonPlt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from shapely import geometry as geom
from shapely.geometry import Polygon, mapping, shape, MultiPolygon
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

import processingMS

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
            #print("mippt", mip_pt)

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
            ax.annotate(f"U:{row['ts_wu']:.0f}, V:{row['ts_wv']:.0f}, E:{row['ts_mipPt']:.2f}", (row['ts_x'], row['ts_y']+4), fontsize=3, color='red')

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
            ax.annotate(f"U:{row['tc_cu']:.0f}, V:{row['tc_cv']:.0f}, E:{row['tc_mipPt']:.2f}", (row['tc_x'], row['tc_y']+4), fontsize=3, color='red')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'TC - Scatter plot of layer {layer}')
        ax.legend()
        plt.savefig(f"sci_plot_of_layer_{layer}.png", dpi=600)
        plt.close()


def plot_baseline(df_baseline_proj, algo, event, particle):
    fig, ax = plt.subplots(figsize=(10, 8))

    # Dictionary to store total mipPt values for each bin
    total_mipPt = {}

    for index, row in df_baseline_proj.iterrows():
        eta_vertices = row['eta_vertices']
        phi_vertices = row['phi_vertices']
        mipPt = row['mipPt']

        # Store total mipPt for the bin
        bin_key = (eta_vertices, phi_vertices)
        total_mipPt[bin_key] = mipPt

    # Normalize total mipPt values for color mapping
    max_mipPt = max(total_mipPt.values())
    norm = Normalize(vmin=0, vmax=max_mipPt)

    # Plot bins with total mipPt values
    for bin_key, mipPt in total_mipPt.items():
        eta_vertices, phi_vertices = bin_key
        poly = PolygonPlt(np.column_stack((eta_vertices, phi_vertices)), closed=True, edgecolor='black')
        ax.add_patch(poly)

        # Set color for the bin based on total mipPt value
        color = ScalarMappable(norm=norm, cmap='viridis').to_rgba(mipPt)
        poly.set_facecolor(color)

        # Add text annotation with total mipPt value
        text_x = np.mean(eta_vertices)
        text_y = np.mean(phi_vertices)
        ax.text(text_x, text_y, f'{mipPt:.1f}', color='black', ha='center', va='center', fontsize=5)

    # Add color bar
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    sm = ScalarMappable(norm=norm, cmap='viridis')
    sm.set_array(list(total_mipPt.values()))
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Total mipPt')

    # Set labels and title
    ax.set_xlabel('Eta')
    ax.set_ylabel('Phi')
    ax.set_title(f'{algo}_{particle}_{event}')

    ax.autoscale()
    plt.savefig(f'{algo}_{particle}_{event}.png', dpi=500)  # Save the plot as an image
    plt.show()

def plot_towers_eta_phi_grid(df_baseline_proj, data_gen, algo, event, particle, subdet):
    fig, ax = plt.subplots(figsize=(10, 8))

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

    # Normalize total mipPt values for color mapping
    max_mipPt = df_baseline_proj['mipPt'].max()
    norm = Normalize(vmin=0, vmax=max_mipPt)

    # Plotting bins with colors and annotations
    for _, row in df_baseline_proj.iterrows():
        eta_vertices = row['eta_vertices']
        phi_vertices = row['phi_vertices']
        mipPt = row['mipPt']

        bin_key = (eta_vertices, phi_vertices)
        color = ScalarMappable(norm=norm, cmap='viridis').to_rgba(mipPt)

        poly = PolygonPlt(np.column_stack((eta_vertices, phi_vertices)), closed=True, edgecolor='black')#closed=True
        ax.add_patch(poly)
        poly.set_facecolor(color)

        text_x = np.mean(eta_vertices)
        text_y = np.mean(phi_vertices)
        ax.text(text_x, text_y, f'{mipPt:.1f}', color='black', ha='center', va='center', fontsize=5)

    #if event != '-1':
        # Overlay gen_eta and gen_phi points
        #print('gen_eta', data_gen['gen_eta'])
        #print('gen_phi', data_gen['gen_phi'])
        #print('gen_eta', data_gen['gen_pt'])
        ax.scatter(data_gen['gen_eta'], data_gen['gen_phi'], marker='x', color='red', label='Gen Points')

    # Add color bar
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    sm = ScalarMappable(norm=norm, cmap='viridis')
    sm.set_array(df_baseline_proj['mipPt'])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Total mipPt')

    # Set labels and title
    ax.set_xlabel('Eta')
    ax.set_ylabel('Phi')
    ax.set_title(f'{algo}_{particle}_{event}')

    ax.autoscale()
    plt.savefig(f'{algo}_{particle}_{event}_{subdet}_SUMMED_NEW.png', dpi=500)  # Save the plot as an image
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

                    # Retrieve mippt value from the hexagon's group
                    mippt = hex_group['ts_mipPt'][()]
                    
                    # Retrieve centroid coordinates from the hexagon's group
                    hex_center_x = hex_group['hex_x_centroid'][()]
                    hex_center_y = hex_group['hex_y_centroid'][()]
                    
                    plt.text(hex_center_x, hex_center_y, f'mipPt: {mippt:.2f}', color='red', ha='center', va='center')
                    
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

def plot_window_with_subwindows(window_bins, eta_min, eta_max, phi_min, phi_max, eta_start, eta_end, phi_start, phi_end):
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot the 12x12 window
    for idx, row in window_bins.iterrows():
        polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='black', facecolor='none')
        ax.add_patch(polygon)

    # Highlight the 3x3 subwindow
    subwindow_bins = window_bins[
        (window_bins['eta_center'] >= eta_start) & (window_bins['eta_center'] < eta_end) &
        (window_bins['phi_center'] >= phi_start) & (window_bins['phi_center'] < phi_end)
    ]

    for idx, row in subwindow_bins.iterrows():
        polygon = plt.Polygon(np.column_stack((row['eta_vertices'], row['phi_vertices'])), edgecolor='red', facecolor='none')
        ax.add_patch(polygon)

    # Display the plot
    plt.xlim(eta_min, eta_max)
    plt.ylim(phi_min, phi_max)
    plt.xlabel('Eta')
    plt.ylabel('Phi')
    plt.title(f'3x3 Subwindow from ({eta_start}, {phi_start}) to ({eta_end}, {phi_end})')
    plt.show()

def plot_eta_phi_resolution(results_df):
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
    plt.title(f'Histogram of Eta Differences\n mu={mu_eta:.2f}, sigma={std_eta:.2f}')
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
    plt.title(f'Histogram of Phi Differences\n mu={mu_phi:.2f}, sigma={std_phi:.2f}')
    plt.savefig('eta_phi_diff_histogram_CEH_event_4875.png')

    plt.tight_layout()
    plt.show()