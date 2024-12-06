# coding: utf-8

_all_ = [ 'fill' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

from bye_splits.utils import common
from data_handle import data_process

import random; random.seed(18)
import numpy as np
import pandas as pd
import h5py
from tqdm import tqdm


def fill(pars, df_gen, df_cl, df_tc, df_ts, **kw):
    """
    Fills split clusters information according to the Stage2 FPGA fixed binning.
    NB: At the moment this Filling step is skipped for Modules Sums treatment  -> root files are read directly
    """

    # Sets output paths for generated clusters and all trigger cells.
    out_gencl = common.fill_path(kw['FillOutGenCl'], **pars)
    outtc_all = common.fill_path(kw['FillOutTcAll'], **pars)
    out_ts = common.fill_path(kw['FillOutTs'], **pars)
    with pd.HDFStore(out_gencl, mode='w') as store_gencl, pd.HDFStore(outtc_all, mode='w') as store_all, pd.HDFStore(out_ts, mode='w') as store_ts:

        #Performs baseline selection on the input data (df_gen and df_cl), storing the results in HDF files.
        df1 = data_process.baseline_selection(df_gen, df_cl, pars['sel'], **kw)
        assert(df1[df1.cl3d_eta<0].shape[0] == 0) #It checks whether there are any clusters (cl3d_eta) with negative values after the baseline selection. If there are any, it raises an AssertionError.
                
        #storeComp['df'] = df1.set_index('event').filter(regex='^[gen.*')
        store_gencl['df'] = df1.filter(regex='^gen.*|event')

        df_3d = df1[:].reset_index()

        df_tc['R'] = np.sqrt(df_tc.tc_x**2 + df_tc.tc_y**2)
        df_tc['Rz'] = df_tc.R / abs(df_tc.tc_z)

        df_ts['R'] = np.sqrt(df_ts.ts_x**2 + df_ts.ts_y**2)
        df_ts['Rz'] = df_ts.R / abs(df_ts.ts_z)
        
        # pandas 'cut' returns np.nan when value lies outside the binning
        
        rzedges = np.linspace(kw['MinROverZ'], kw['MaxROverZ'], num=kw['NbinsRz']+1)
        phiedges = np.linspace(kw['MinPhi'], kw['MaxPhi'], num=kw['NbinsPhi']+1)
        phiedgesTs = np.linspace(kw['MinPhi'], kw['MaxPhi'], num=kw['NbinsPhiTs']+1)
        etaedges = np.linspace(kw['EtaMin'], kw['EtaMax'], num=kw['NbinsEta']+1)

        df_tc['Rz_bin'] = pd.cut(df_tc.Rz, bins=rzedges, labels=False)
        df_tc['tc_phi_bin'] = pd.cut(df_tc.tc_phi, bins=phiedges, labels=False)
         #df_ts['ts_Rz_bin'] = pd.cut(df_ts.Rz, bins=rzedges, labels=False)
        df_ts['ts_phi_bin'] =pd.cut(df_ts.ts_phi, bins=phiedgesTs, labels=False)
        df_ts['ts_eta_bin'] =pd.cut(df_ts.ts_eta, bins=etaedges, labels=False)


        # Check for NaN values in 'Rz' and 'phi' columns for tc
        nan_check_rz = df_tc['Rz_bin'].isna().any()
        nan_check_phi = df_tc['tc_phi_bin'].isna().any()
        #print("NaN check for 'Rz' column:", nan_check_rz)
        #print("NaN check for 'tc_phi' column:", nan_check_phi)

        # Check for NaN values in 'eta' and 'phi' columns for ts
        nan_check_etaTS = df_ts['ts_eta_bin'].isna().any()
        nan_check_phTS = df_ts['ts_phi_bin'].isna().any()
        #print("NaN check for 'ts_eta ' column:", nan_check_etaTS)
        #print("NaN check for 'ts_phi' column:", nan_check_phTS)
        # Print the rows with NaN values in specific columns
        nan_rows = df_ts[df_ts['ts_eta_bin'].isna()]
        #print("Rows with NaN values in 'Rz' or 'tc_phi' columns:")
        #print(nan_rows)
        

        nansel = (pd.isna(df_tc.Rz_bin)) & (pd.isna(df_tc.tc_phi_bin))
        nanselts = (pd.isna(df_ts.ts_phi_bin)) | (pd.isna(df_ts.ts_eta_bin))
        df_tc = df_tc[~nansel]
        df_ts = df_ts[~nanselts]

       # data_process.plot_2d(df_tc_event, 'tc_phi_bin' , 'Rz_bin', phiedges, rzedges,'Phi','Rz', '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/tasks/tc_plot_event_phi_Rz.png')
       # data_process.plot_2d(df_ts_event, 'ts_phi_bin' , 'ts_eta_bin', phiedges, etaedges,'Phi', 'Eta', '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/tasks/tsum_plot_event_phi_eta.png')
       # data_process.plot_2d(df_ts_event, 'ts_phi_bin' , 'ts_Rz_bin', phiedges, etaedges,'Phi', 'Rz', '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/tasks/tsum_plot_event_phi_Rz.png')

        #data_process.plot_2d(df_tc, 'tc_phi_bin' , 'Rz_bin', phiedges, rzedges,'Phi','Rz', '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/tasks/tc_plot_phi_Rz.png')
        #data_process.plot_2d(df_ts, 'ts_phi_bin' , 'ts_eta_bin', phiedgesTs, etaedges,'Phi', 'Eta', '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/tasks/tsum_plot_phi_eta_binNew.png')
       # data_process.plot_2d(df_ts, 'ts_phi_bin' , 'ts_Rz_bin', phiedges, etaedges,'Phi', 'Rz', '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/tasks/tsum_plot_phi_Rz.png')

        store_all[kw['FesAlgo'] + '_3d'] = df_3d
        store_ts['ts']=df_ts

        """
        hexagon = [
        [1, 0, 3],
        [0.5, np.sqrt(3)/2, 3],
        [-0.5, np.sqrt(3)/2, 3],
        [-1, 0, 0],
        [-0.5, -np.sqrt(3)/2, 3],
        [0.5, -np.sqrt(3)/2, 3],
        ]

        hexagon_transformed = []

        # Convert coordinates to cylindrical coordinates (phi, eta)
        for vertex in hexagon:
            phi, eta = data_process.cartesian_to_cylindrical(*vertex)
            hexagon_transformed.append([phi, eta])

        # Plotting
        data_process.plot_hexagons(hexagon, hexagon_transformed)
        """
        data_process.plot_grids_comparison()
        dfs = (df_3d, df_tc, df_ts) 
            
    ### Event Processing ######################################################
    outfill  = common.fill_path(kw['FillOut'], **pars)
    with h5py.File(outfill, mode='w') as store, pd.HDFStore(outtc_all, mode='w') as store_all:
        group_tot = None
        df_3d, df_tc, df_ts = dfs

        for ev in tqdm(df_tc['event'].unique().astype('int')):
            ev_tc = df_tc[df_tc.event == ev]
            ev_3d = df_3d[df_3d.event == ev]
            ev_ts = df_ts[df_ts.event == ev]
            if ev_3d.empty or ev_tc.empty or ev_ts.empty :
                continue
            
            keep_tc = ['tc_phi_bin', 'Rz_bin', 'tc_layer', 'tc_mipPt', 'tc_pt',
                       'tc_wu', 'tc_wv', 'tc_cu',  'tc_cv', 'tc_x', 'tc_y', 'tc_z',
                       'tc_eta', 'tc_phi', 'gen_eta', 'gen_phi', 'tc_subdet']
            ev_tc = ev_tc.filter(items=keep_tc)
            wght_f = lambda pos: ev_tc.tc_mipPt*pos/np.abs(ev_tc.tc_z)
            ev_tc['wght_x'] = wght_f(ev_tc.tc_x)
            ev_tc['wght_y'] = wght_f(ev_tc.tc_y)
  
            with common.SupressSettingWithCopyWarning():
                ev_3d['cl3d_rz'] = common.calcRzFromEta(ev_3d.loc[:,'cl3d_eta'])
                ev_3d['gen_rz']  = common.calcRzFromEta(ev_3d.loc[:,'gen_eta'])

            cl3d_rz  = ev_3d['cl3d_rz'].unique()
            cl3d_phi = ev_3d['cl3d_phi'].unique()
            cl3d_eta = ev_3d['cl3d_eta'].unique()
            cl3d_en  = ev_3d['cl3d_en'].unique()

            store_str = kw['FesAlgo'] + '_' + str(ev) + '_cl'
            cl3d_info = {'cl3d_eta': cl3d_eta, 'cl3d_phi': cl3d_phi,
                         'cl3d_rz': cl3d_rz, 'cl3d_en': cl3d_en}
            store[store_str] = list(cl3d_info.values())
            store[store_str].attrs['columns'] = list(cl3d_info.keys())
            store[store_str].attrs['doc'] = 'CMSSW cluster info.'
            
            gen_rz = ev_3d['gen_rz'].unique()
            gen_phi = ev_3d['gen_phi'].unique()
            ev_3d = ev_3d.drop(['cl3d_rz', 'cl3d_eta', 'cl3d_phi'], axis=1)
            if len(gen_rz) != 1 or len(gen_phi) != 1:
                mes = 'Impossible! Rz: {} | Phi: {}'.format(gen_rz, gen_phi)
                raise RuntimeError(mes)
            
            group = ev_tc.groupby(['Rz_bin', 'tc_phi_bin'], as_index=False)
            cols_keep = ['Rz_bin', 'tc_phi_bin', 'tc_mipPt', 'wght_x', 'wght_y']
  
            group = group.sum()[cols_keep]
            group.loc[:, ['wght_x', 'wght_y']] = group.loc[:, ['wght_x', 'wght_y']].div(group.tc_mipPt, axis=0)

            store_str = kw['FesAlgo'] + '_' + str(ev) + '_ev'
            store[store_str] = group.to_numpy()
            store[store_str].attrs['columns'] = cols_keep
            store[store_str].attrs['doc'] = 'R/z vs. Phi histo Info'

            cols_keep = ['Rz_bin', 'tc_phi_bin', 'tc_layer', 'tc_mipPt', 
                         'tc_pt', 'tc_wu', 'tc_wv', 'tc_cu',  'tc_cv',
                         'tc_x',  'tc_y',  'tc_z',  'tc_eta', 'tc_phi', 'tc_subdet']
            ev_tc = ev_tc[cols_keep]
            if ev_tc.empty:
                continue

            store_str = kw['FesAlgo'] + '_' + str(ev) + '_tc'
            store[store_str] = ev_tc.to_numpy()
            store[store_str].attrs['columns'] = cols_keep
            store[store_str].attrs['doc'] = 'Trigger Cells Info'
            
            #NEW TS
            cols_keep = ['ts_phi_bin', 'ts_eta_bin', 'ts_layer', 'ts_mipPt', 
                         'ts_pt', 'ts_wu', 'ts_wv','ts_x',  'ts_y',  'ts_z',
                         'ts_eta', 'ts_phi', 'ts_energy', 'ts_subdet']
            ev_ts = ev_ts[cols_keep]
            if ev_ts.empty:
                continue

            store_str = kw['FesAlgo'] + '_' + str(ev) + '_ts'
            store[store_str] = ev_ts.to_numpy()
            store[store_str].attrs['columns'] = cols_keep
            store[store_str].attrs['doc'] = 'Module Sums Info'

            ########################

            ev_tc_all = df_tc[df_tc.event == ev]
            ev_tc_all = ev_tc_all.reset_index().drop(['entry', 'subentry', 'event'], axis=1)
            all_keep = ['tc_layer', 'tc_mipPt', 'tc_pt', 'tc_energy',
                        'tc_x', 'tc_y', 'tc_z', 'tc_eta', 'tc_phi','tc_subdet']
            all_tcs = ev_tc_all.filter(items=all_keep)
            if ev_tc_all.empty:
                continue
            store_all[store_str] = all_tcs

            #print("\nall tcs:")
            #print(all_tcs)
            #print("Length of all_tcs:", len(all_tcs))

            if group_tot is not None:
                group_tot = group[:]
            else:
                group_tot = pd.concat((group_tot, group[:]), axis=0)

        return group_tot

if __name__ == "__main__":
    import argparse
    from bye_splits.utils import params, parsing

    parser = argparse.ArgumentParser(description='Filling standalone step.')
    parsing.add_parameters(parser)
    FLAGS = parser.parse_args()
    assert (FLAGS.sel in ('splits_only', 'no_splits', 'all') or
            FLAGS.sel.startswith('above_eta_'))
    df_gen, df_cl, df_tc, df_ts = data_process.get_data_reco_chain_start_ModSums(nevents=-1, reprocess=True, tag='chain', particles='pions', pu=0, event=None)
    fill_d = params.read_task_params('fill')
    fill(vars(FLAGS), df_gen, df_cl, df_tc, df_ts, **fill_d)