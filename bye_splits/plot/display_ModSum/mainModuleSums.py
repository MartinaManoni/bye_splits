# coding: utf-8
# python3 mainModuleSums.py --event -1 --geom V16 --algo baseline --particle pions --subdet 5

_all_ = [ ]

import os
import sys

parent_dir = os.path.abspath(__file__ + 4 * '/..')
sys.path.insert(0, parent_dir)

import argparse
import processingMS
import resolutionMS
import helperMS
import geometryMS
import json
import warnings
warnings.filterwarnings("ignore", category=UserWarning, message=".*subnormal.*")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Interactive Grid Comparison")

    parser.add_argument("--subdet", type=int, default=1, help="1: CEE (has only silicon layers), 2: CEH - only silicon part, 3: CEH - only scint part, 4: CEH, all layers, 5: CEE + CEH")
    parser.add_argument("--event", default='5492', help="Select event number or -1 for all events")
    parser.add_argument("--n", type=int, default=None, help="Process n events (random ordering)")
    parser.add_argument("--algo", default='8towers', help="Select algorithm (baseline, area_overlap, 8towers, 16towers)")
    parser.add_argument("--particle", default='photons', help="Select particle type (photons, pions or neutrinos)")
    parser.add_argument("--geom", default='V16', help="Select the CMSSW geometry (V11 or V16)")
    return parser.parse_args()


def main(subdet, event, particle, algo, n, geom):
    process = processingMS.Processing() 
    #algorithms = algosMS.Algorithms()
    #resolution = resolutionMS.Resolution()
    #helper = helperMS.Helper()
    geometry = geometryMS.GeometryData()
    print("entering main")

    save_tower_bins = False

    #file_path = f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/DoublePhotonsPU0_3k_V11/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
    if geom=='V11':
        print("using V11 geom")
        file_path = (
        '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/'
        'DoublePhotonsPU0_hadd_123_energy/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
        )
    elif geom=='V16' and particle == 'photons':
        file_path = (
        '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/'
        'SinglePhotonPU0V16/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
        )
        root_file =('/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples/SinglePhotonPU0V16.root')

    elif geom=='V16' and particle == 'pions':
        file_path = (
        '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/'
        'SinglePionPU0V16/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
        )
        root_file =('/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples/SinglePionPU0V16.root')

    elif geom=='V16' and particle == 'neutrinos':
        root_file =('/data_CMS/cms/manoniL1HGCAL/ntupleV16Production/MinBias_TuneCP5_14TeV-pythia8_Phase2Fall22DRMiniAOD-PU200_125X_mcRun4_realistic_v2-v1/skimmed_ntuples/Ntuple_1.root')


    # Skip get_gen_particles for neutrinos
    if particle == 'neutrinos':
        print("Skipping get_gen_particles for neutrinos")
        df_gen, events = None, None  # Set these to None since gen variables are unavailable
    else:
        df_gen, events = process.get_gen_particles(root_file, n, event)
        print("df_gen", df_gen)
        print("events", events)

    # Create the specific data frame for the particle
    if events is None or len(events) == 0:
        selected_events = None
    else:
        selected_events = events

    df_specific = process.read_root_and_create_dataframe(
    root_file, subdet, selected_events= selected_events)


    #df_specific = process.read_root_and_create_dataframe(root_file, subdet, events)
    #print("df_specific", df_specific)
    #print("df_specific", df_specific.columns)

    #helper.read_hdf5_structure(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/photons_manoni/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5')
    #helper.read_all_block0_values(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/photons_manoni/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5')

    bin_geojson_filename = '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_with_arcs.geojson'
    hex_geojson_filename = '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/hexagons_CMSSW.geojson'
    scint_geojson_filename = '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/scint_modules_geo.geojson'

    #overlap = process.eval_hex_bin_overlap(data, bin_geojson_filename,  hdf5_filename)

    #hexagon_info_df = process.eval_hex_bin_overlap_OK(data, bin_geojson_filename)

    #cProfile.run('process.eval_hex_bin_overlap(data, bin_geojson_filename,  hdf5_filename)')

    initial_kw = {
        'NbinsEta': 20,
        'NbinsPhi': 72,
        'MinPhi': -3.14159,
        'MaxPhi': +3.14159,
        'EtaMin': 1.305,
        'EtaMax': 3.045
    }


    #bins_data, hexagons_data, scint_data = geometry.read_geojson_files(bin_geojson_filename, hex_geojson_filename, scint_geojson_filename)
    #plotMS.plot_full_geom(bins_data, hexagons_data, scint_data, 'plot_layers', plot_type='all')

    if save_tower_bins:
        process.create_and_save_tower_bins(initial_kw, geom) #create and save tower bins

    process.ModSumToTowers(initial_kw, df_specific , subdet, event, particle, algo, bin_geojson_filename, hex_geojson_filename, df_gen, geom)

    #geometry.save_bin_geo(towers_bins, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_with_arcs.geojson', f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_only_vertices.geojson')
    #geometry.save_bin_hex(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/hexagons_CMSSW.geojson')
    #geometry.save_scint_mod_geo(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/scint_modules_geo.geojson')

if __name__ == '__main__':
    args = parse_arguments()
    main(args.subdet, args.event, args.particle, args.algo, args.n, args.geom)

    

