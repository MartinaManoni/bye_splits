# coding: utf-8

_all_ = ['baseline_selection', 'EventDataParticle']

import os
import sys

parent_dir = os.path.abspath(__file__ + 2 * "/..")
sys.path.insert(0, parent_dir)

import numpy as np
import pandas as pd
import yaml
import h5py
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColorBar, BasicTicker, LinearColorMapper
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

import bye_splits
from bye_splits.utils import common
from utils import params
from data_handle.geometry import GeometryData
from data_handle.event import EventData
from data_handle.data_input import InputData

def baseline_selection(df_gen, df_cl, sel, **kw):
    data = pd.merge(left=df_gen, right=df_cl, how='inner', on='event')
    nin = data.shape[0]
    data = data[(data.gen_eta>kw['EtaMin']) & (data.gen_eta<kw['EtaMax'])] #Filters the data based on the range of gen_eta

    #if sel starts with 'above_eta_', it further filters the data to include only events with gen_eta greater than the specified threshold.
    if sel.startswith('above_eta_'):
        data = data[data.gen_eta > float(sel.split('above_eta_')[1])]
        return data
    
    #Calculates the energy resolution (enres) by subtracting gen_en from cl3d_en and normalizing by gen_en.
    with common.SupressSettingWithCopyWarning():
        data['enres'] = data.cl3d_en - data.gen_en #enres = energy resolution
        data.enres /= data.gen_en

    nansel = pd.isna(data['enres'])
    nandf = data[nansel]
    nandf['enres'] = 1.1
    data = data[~nansel]
    data = pd.concat([data,nandf], sort=False)

        #If sel is 'splits_only', it selects events with clusters having enres (energy resolution) less than kw['EnResSplits'].
    if sel == 'splits_only':
        # select events with splitted clusters (enres < energy cut)
        # if an event has at least one cluster satisfying the enres condition,
        # all of its clusters are kept (this eases comparison with CMSSW)
        evgrp = data.groupby(['event'], sort=False)
        multiplicity = evgrp.size()
        bad_res = (evgrp.apply(lambda grp: np.any(grp['enres'] < kw['EnResSplits']))).values
        bad_res_mask = np.repeat(bad_res, multiplicity.values)
        data = data[bad_res_mask]

    elif sel == 'no_splits':
        data = data[(data.gen_eta > kw['EtaMinStrict']) &
                    (data.gen_eta < kw['EtaMaxStrict'])]
        evgrp = data.groupby(['event'], sort=False)
        multiplicity = evgrp.size()
        good_res = (evgrp.apply(lambda grp: np.all(grp['enres'] > kw['EnResNoSplits']))).values
        good_res_mask = np.repeat(good_res, multiplicity.values)
        data = data[good_res_mask]

    #If sel is 'all', it does not perform any additional selection.
    elif sel == 'all':
        pass

    else:
        m = 'Selection {} is not supported.'.format(sel)
        raise ValueError(m)

    nout = data.shape[0]
    eff = (nout / nin) * 100
    print("The baseline selection has a {}% efficiency: {}/{}".format(np.round(eff,2), nout, nin))
    return data   #Returns the final DataFrame after applying all the specified selection criteria


#The get_data_reco_chain_start function orchestrates the retrieval of
# event data for generated particles, clusters, and trigger cells.
# It uses the EventDataParticle class to handle the specifics of
# loading data for different particle types based on the provided parameters.
# The relevant columns are selected and returned as separate datasets
def get_data_reco_chain_start(nevents=500, reprocess=False, tag='chain', particles='photons', pu=0, event=None): #qui event è  None quindi di default vengono slezionati random events
    """Access event data."""

    #Creates an instance of the EventDataParticle class, specifying the particle type (particles), pileup (pu), tag, and whether to reprocess the data (reprocess)
    data_particle = EventDataParticle(particles, pu, tag, reprocess)
    if event is None: #if event is None, random events are retrieved using the provide_random_events method of the EventDataParticle instance.
        ds_all, events = data_particle.provide_random_events(n=nevents)
        # ds_all = data_particle.provide_events(events=[170004, 170015, 170017, 170014])
    else: #if event is specified, data for that specific event is retrieved using the provide_event method.
        ds_all = data_particle.provide_event(event, merge=False)
        events = event

    if ds_all["gen"].empty:
        raise RuntimeError("No events in the parquet file.")

    tc_keep = {
        "event": "event",
        "good_tc_waferu": "tc_wu",
        "good_tc_waferv": "tc_wv",
        "good_tc_cellu" : "tc_cu",
        "good_tc_cellv" : "tc_cv",
        "good_tc_layer" : "tc_layer",
        "good_tc_pt"    : "tc_pt",
        "good_tc_mipPt" : "tc_mipPt",
        "good_tc_energy": "tc_energy",
        "good_tc_x"     : "tc_x",
        "good_tc_y"     : "tc_y",
        "good_tc_z"     : "tc_z",
        "good_tc_eta"   : "tc_eta",
        "good_tc_phi"   : "tc_phi",
    }

    ds_tc = ds_all["tc"]
    ds_tc = ds_tc[tc_keep.keys()]
    ds_tc = ds_tc.rename(columns=tc_keep)

    gen_keep = {
        "event": "event",
        "good_genpart_exeta" : "gen_eta",
        "good_genpart_exphi" : "gen_phi",
        "good_genpart_energy": "gen_en",
        "good_genpart_pt"    : "gen_pt",
    }
    ds_gen = ds_all["gen"]
    ds_gen = ds_gen.rename(columns=gen_keep)

    cl_keep = {
        "event": "event",
        "good_cl3d_eta"   : "cl3d_eta",
        "good_cl3d_phi"   : "cl3d_phi",
        "good_cl3d_id"    : "cl3d_id",
        "good_cl3d_energy": "cl3d_en",
        "good_cl3d_pt"    : "cl3d_pt",
    }
    ds_cl = ds_all["cl"]
    ds_cl = ds_cl.rename(columns=cl_keep)
    return ds_gen, ds_cl, ds_tc
            
def EventDataParticle(particles, pu, tag, reprocess, logger=None):
    """Factory for EventData instances of different particle types"""
    with open(params.CfgPath, "r") as afile:
        cfg = yaml.safe_load(afile)
        if particles is None:
            particles = cfg["selection"]["particles"]
        if particles not in ("photons", "electrons", "pions"):
            raise ValueError("{} are not supported.".format(particles))
        defevents = cfg["defaultEvents"][f"PU{pu}"][particles]
        #print("print 3")
        indata = InputData()
        #print("print 4")
        indata.path = cfg["io"][f"PU{pu}"][particles]["file"]
        indata.adir = cfg["io"][f"PU{pu}"][particles]["dir"]
        indata.tree = cfg["io"][f"PU{pu}"][particles]["tree"]
        #print("print 5")
    tag = particles + "_" + f"PU{pu}" + "_" + tag 
    #print("print 6")
    return EventData(indata, tag, defevents, reprocess, logger)

def get_data_reco_chain_start_ModSums(nevents=500, reprocess=False, tag='chain', particles='photons', pu=0, event=None): #qui event è  None quindi di default vengono slezionati random events
    """Access event data."""
    print("entered get data reco")
    #Creates an instance of the EventDataParticle class, specifying the particle type (particles), pileup (pu), tag, and whether to reprocess the data (reprocess)
    data_particle = EventDataParticle(particles, pu, tag, reprocess)
    if event is None: #if event is None, random events are retrieved using the provide_random_events method of the EventDataParticle instance.
        ds_all, event = data_particle.provide_random_events(n=nevents)
        # ds_all = data_particle.provide_events(events=[170004, 170015, 170017, 170014])
    else: #if event is specified, data for that specific event is retrieved using the provide_event method.
        ds_all = data_particle.provide_event(event, merge=False)
        events = event

    if ds_all["gen"].empty:
        raise RuntimeError("No events in the parquet file.")

    tc_keep = {
        "event": "event",
        "good_tc_waferu": "tc_wu",
        "good_tc_waferv": "tc_wv",
        "good_tc_cellu" : "tc_cu",
        "good_tc_cellv" : "tc_cv",
        "good_tc_layer" : "tc_layer",
        "good_tc_pt"    : "tc_pt",
        "good_tc_mipPt" : "tc_mipPt",
        "good_tc_energy": "tc_energy",
        "good_tc_x"     : "tc_x",
        "good_tc_y"     : "tc_y",
        "good_tc_z"     : "tc_z",
        "good_tc_eta"   : "tc_eta",
        "good_tc_phi"   : "tc_phi",
        "good_tc_subdet": "tc_subdet",
    }

    ds_tc = ds_all["tc"]
    ds_tc = ds_tc[tc_keep.keys()]
    ds_tc = ds_tc.rename(columns=tc_keep)

    ts_keep = {
        "event": "event",
        "good_ts_waferu": "ts_wu",
        "good_ts_waferv": "ts_wv",
        "good_ts_layer" : "ts_layer",
        "good_ts_pt"    : "ts_pt",
        "good_ts_mipPt" : "ts_mipPt",
        "good_ts_energy": "ts_energy",
        "good_ts_x"     : "ts_x",
        "good_ts_y"     : "ts_y",
        "good_ts_z"     : "ts_z",
        "good_ts_eta"   : "ts_eta",
        "good_ts_phi"   : "ts_phi",
        "good_ts_subdet": "ts_subdet",

    }

    ds_ts = ds_all["tsum"]
    ds_ts = ds_ts[ts_keep.keys()]
    ds_ts = ds_ts.rename(columns=ts_keep)

    gen_keep = {
        "event": "event",
        "good_genpart_exeta" : "gen_eta",
        "good_genpart_exphi" : "gen_phi",
        "good_genpart_energy": "gen_en",
        "good_genpart_pt"    : "gen_pt",
    }
    ds_gen = ds_all["gen"]
    ds_gen = ds_gen.rename(columns=gen_keep)

    #print("EVENT 5492", ds_gen[ds_gen.event == 5492])

    cl_keep = {
        "event": "event",
        "good_cl3d_eta"   : "cl3d_eta",
        "good_cl3d_phi"   : "cl3d_phi",
        "good_cl3d_id"    : "cl3d_id",
        "good_cl3d_energy": "cl3d_en",
        "good_cl3d_pt"    : "cl3d_pt",
    }
    ds_cl = ds_all["cl"]
    ds_cl = ds_cl.rename(columns=cl_keep)



    return ds_gen, ds_cl, ds_tc, ds_ts


def plot_2d(df, x_col, y_col, x_edges, y_edges, xlabel, ylabel, outname):
    """
    Plot a grid in phi and rz and visualize 2D plots with events in each bin.

    Parameters:
    - df: DataFrame containing event data.
    - phi_col: Column name for the phi values.
    - rz_col: Column name for the rz values.
    - phi_edges: Bin edges for phi.
    - rz_edges: Bin edges for rz.
    """

    # Create a 2D histogram
    #hist, x_edges, y_edges = np.histogram2d(df[phi_col], df[rz_col], bins=[phi_edges, rz_edges])
    hist, _, _ = np.histogram2d(df[x_col], df[y_col], bins=[len(x_edges) - 1, len(y_edges) - 1])

    # Create phi and rz grid
    x_grid, y_grid = np.meshgrid(x_edges, y_edges)

    # Plot the 2D histogram
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(x_grid, y_grid, hist.T, cmap='viridis', shading='auto')
    plt.colorbar(label='Number of Events')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('2D Plot')
    plt.grid(True)
    plt.savefig(outname)