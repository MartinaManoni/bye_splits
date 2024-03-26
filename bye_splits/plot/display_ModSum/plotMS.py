# coding: utf-8

_all_ = []

import os
from pathlib import Path
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)

from bye_splits.plot.display_plotly import yaml, np, pd, go, dcc

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from scipy.spatial.distance import cdist

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

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
    fig.add_trace(go.Scatter(
        x=df_bin['centroid_X'],
        y=df_bin['centroid_Y'],
        mode='markers',
        marker=dict(color='black', size=1),
        name='Center Points'
    ))    

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


def plot_single_event(grid_params, hexagon_data, df_bin, output_path):
    
    sample_df = processingMS.Processing().create_grid_df(grid_params)
    test_fig = plot_grid(sample_df, df_bin, grid_params)

    grouped_data = hexagon_data.groupby('ts_layer')

    # Iterate through the groups
    for ts_layer_value, group in grouped_data:
        # Iterate through the rows of the group
        for index, row in group.iterrows():
            hex_x = row['hex_x']
            hex_y = row['hex_y']

            # Close the hexagon by repeating the first vertex at the end
            hex_x.append(hex_x[0])
            hex_y.append(hex_y[0])

            # Add the hexagon trace to the existing Plotly figure
            #test_fig.add_trace(px.line(x=hex_x, y=hex_y, line_shape='linear').data[0])
              # Add the hexagon trace to the existing Plotly figure
            test_fig.add_trace(go.Scatter(
                x=hex_x,
                y=hex_y,
                mode='lines',
                #fill='toself',
                #fillcolor='red'
                #line=dict(color='red'),
                #marker=dict(size=12, color=row['ts_mipPt'], colorscale='Viridis', colorbar=dict(title='mipPt'))
                #marker=dict(size=12,color=[row['ts_mipPt']] * len(hex_x),colorscale='Viridis', colorbar=dict(title='mipPt'))
            ))
    test_fig.write_image(output_path)
    return test_fig 

def plot_single_event_proj(grid_params, hexagon_proj, df_bin, output_path2):
    grid_df = processingMS.Processing().create_grid_df(grid_params)
    test_fig = plot_grid(grid_df, df_bin, grid_params)

    # Iterate through the rows of aggregated_df
    for index, row in hexagon_proj.iterrows():
        hex_x = row['hex_x']
        hex_y = row['hex_y']
        mip_pt = row['ts_mipPt']
        ts_wu = row['ts_wu']
        ts_wv = row['ts_wv']
        #print("mippt", mip_pt)

        # Define fill color based on mip_pt value
        fill_color = 'rgba(0,0,255,0.2)'  # Default fill color
        '''if mip_pt is not None:
            # Adjust the fill color based on the mip_pt value
            fill_color = f'rgba(0,0,255,{mip_pt/max(hexagon_proj["ts_mipPt"])})'
            print("mippt MAX", max(hexagon_proj["ts_mipPt"]))
            print("mippt", fill_color)'''
        
        if mip_pt is not None and mip_pt / max(hexagon_proj["ts_mipPt"]) > 0.0001:
            # Adjust the fill color based on the mip_pt value
            max_mip_pt = max(hexagon_proj["ts_mipPt"])
            alpha = mip_pt / max_mip_pt
            fill_color = f'rgba(0, 0, 255, {alpha})'
            #print("mippt MAX:", max_mip_pt)
            #print("mippt:", fill_color)
        else:
            fill_color = 'rgba(0, 0, 255, 0)'

        hex_x_closed = hex_x + (hex_x[0],)
        hex_y_closed = hex_y + (hex_y[0],)

        # Add the hexagon trace to the existing Plotly figure
        test_fig.add_trace(go.Scatter(
            x=hex_x_closed,
            y=hex_y_closed,
            mode='lines',
            line=dict(color='red'),  # Adjust the color as needed
            fill='toself',
            fillcolor=fill_color,
            marker=dict(size=12, color=mip_pt, colorscale='Viridis', colorbar=dict(title='mipPt')),
            hoverinfo='text',
            text=f'mipPt: {mip_pt}'
        ))

        # Add the hexagon trace to the existing Plotly figure
        test_fig.add_trace(go.Scatter(
            x=hexagon_proj['ts_x'],
            y=hexagon_proj['ts_y'],
            mode='markers',
            marker=dict(color='red', size=4),  # Adjust the color as needed
            name=f'hex_centroid'
        ))

        test_fig.add_trace(go.Scatter(
            x=hexagon_proj['ts_x'],
            y=hexagon_proj['ts_y'],
            mode='text',
            text=[f'{wu}, {wv}' for wu, wv in zip(hexagon_proj['ts_wu'], hexagon_proj['ts_wv'])],
            textposition='middle center',  # Set the position of the text to the middle of each point
            textfont=dict(size=7, color='black'),  # Adjust the font size and color as needed
        ))

        '''for tsx, tsy, wu, wv in zip(hexagon_proj['ts_x'], hexagon_proj['ts_y'], hexagon_proj['ts_wu'], hexagon_proj['ts_wv']):
            test_fig.add_annotation(
                x=tsx,
                y=tsy,
                text=f'{wu}, {wv}',
                showarrow=False,
                font=dict(size=10, color='black'),  # Adjust the font size and color as needed
                xshift=0,  # Adjust the x-shift if needed
                yshift=0   # Adjust the y-shift if needed
            )'''


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


def plot_grid_area_overlap(merged_df, algo, subdet, event, particle):
    fig, ax = plt.subplots(figsize=(10, 8))

    # Dictionary to store cumulative mipPt values for each bin
    cumulative_mipPt = {}

    for index, row in merged_df.iterrows():
        for bin_info in row['bins_overlapping']:
            eta_vertices = tuple(bin_info['eta_vertices'])
            phi_vertices = tuple(bin_info['phi_vertices'])
            ts_mipPt = bin_info['mipPt']

            # Update cumulative mipPt for the bin
            bin_key = (eta_vertices, phi_vertices)   #here change key?
            if bin_key not in cumulative_mipPt:
                cumulative_mipPt[bin_key] = 0
            cumulative_mipPt[bin_key] += ts_mipPt
            '''if cumulative_mipPt[bin_key] !=0:
                print("bin key", bin_key)    
                print("here", cumulative_mipPt[bin_key])'''

    # Normalize cumulative mipPt values for color mapping
    max_cumulative_mipPt = max(cumulative_mipPt.values())
    norm = Normalize(vmin=0, vmax=max_cumulative_mipPt)

    # Plot bins with cumulative mipPt values
    for bin_key, cumulative_mip in cumulative_mipPt.items():
        eta_vertices, phi_vertices = bin_key
        poly = Polygon(np.column_stack((eta_vertices, phi_vertices)), closed=True, edgecolor='black')
        ax.add_patch(poly)

        # Set color for the bin based on cumulative mipPt value
        color = ScalarMappable(norm=norm, cmap='viridis').to_rgba(cumulative_mip)
        poly.set_facecolor(color)

        # Add text annotation with cumulative mipPt value
        text_x = np.mean(eta_vertices)
        text_y = np.mean(phi_vertices)
        ax.text(text_x, text_y, f'{cumulative_mip:.1f}', color='black', ha='center', va='center', fontsize=5)

    # Add color bar
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    sm = ScalarMappable(norm=norm, cmap='viridis')
    sm.set_array(list(cumulative_mipPt.values()))
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('mipPt')

    # Set labels and title
    ax.set_xlabel('Eta')
    ax.set_ylabel('Phi')
    ax.set_title(f'{algo}_{subdet}_{particle}_{event}')

    ax.autoscale()
    plt.savefig(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/{algo}_{subdet}_{particle}_{event}.png', dpi=500)


def plot_shifted_modules(silicon_df_proc, shifted_df, df_geom, df_ts, layer_number):
    '''function to plot the centers of module sums for byesplit geometry 
    and CMSSSW, the two geometry have a relative shift'''

    before_shift = silicon_df_proc[(silicon_df_proc['ts_layer'] == layer_number) & (silicon_df_proc['ts_eta']> 0) ]
    after_shift = shifted_df[(shifted_df['ts_layer'] == layer_number)]   
    df_geom_m = df_geom[(df_geom['ts_layer']==layer_number)]
    df_ts = df_ts[(df_ts['ts_layer']==layer_number) & (df_ts['ts_eta']> 0) ]

    plt.figure(figsize=(8, 6))
    plt.scatter(before_shift['wx_center'], before_shift['wy_center'], label='Byesplit geometry (after merge)', color='red', s=1)
    plt.scatter(before_shift['ts_x'], before_shift['ts_y'], label='ntuples data (after merge)', color='black', s=1)
    #plt.scatter(df_ts['ts_x'], df_ts['ts_y'], label='ntuples data (before merge)', color='cyan', s=1)
    #plt.scatter(df_geom_m['wx_center'], df_geom_m['wy_center'], label='CMSSW geom marco (before emerge)', color='green', s=0.5)

    # Plot wafer U and wafer V for wx_center/wy_center (byesplit geometry)
    for idx, row in before_shift.iterrows():
        plt.text(row['wx_center'], row['wy_center']+2, f"U:{row['ts_wu']}, V:{row['ts_wv']}", fontsize=3, color='red')

    # Plot wafer U and wafer V for ts_x/ts_y  (ntuples data after data/geometry merge)
    for idx, row in before_shift.iterrows():
        plt.text(row['ts_x'], row['ts_y'], f"U={row['ts_wu']}, V={row['ts_wv']}", fontsize=3)    

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
    plt.savefig(f"/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/Layer{layer_number}_OKAY_hex.png", dpi = 400)

def plot_group_positions(final_projected_df):
    # Define colors for each group
    group_colors = {
        'CEE': 'blue',
        'CEH_odd': 'green',
        'CEH_even': 'red'
        #'CEH_odd_2': 'orange',
       #'CEH_even_2': 'cyan'
    }
    # Create a subplot for each group
    for group, color in group_colors.items():
        group_df = final_projected_df[final_projected_df['group'] == group]
        mask = (group_df["ts_wu"] == 0) & (group_df["ts_wv"] == 6)
        plt.figure(figsize=(8, 6))
        
        for hex_x, hex_y in zip(group_df['hex_x'], group_df['hex_y']):
            plt.plot(hex_x + (hex_x[0],), hex_y + (hex_y[0],), color=color, alpha=0.5)  # Close the loop

        #plt.scatter(group_df['ts_x'], group_df['ts_y'], color=color, s=1)
        plt.scatter(group_df['wx_center'], group_df['wy_center'], color=color, s=1)
        plt.title(f"Position of Entries for {group} with Annotations")
        plt.xlabel('ts_x')
        plt.ylabel('ts_y')
        
        # Annotate each point with ts_wu and ts_wv values
        for index, row in group_df.iterrows():
            plt.annotate(f"U:{row['ts_wu']}, V:{row['ts_wv']}", (row['wx_center'], row['wy_center']), fontsize=3)
            plt.annotate(f"E:{row['ts_mipPt']:.2f}", (row['wx_center'], row['wy_center']-3), fontsize=3, color=color)
        
        print(f"DF GROUP {group}",group_df)
        print(f"lenght {group}", len(group_df))
        # Save the plot as PNG
        plt.savefig(f"{group}_positions_pion_NEW.png", dpi=400)
        plt.close()


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
