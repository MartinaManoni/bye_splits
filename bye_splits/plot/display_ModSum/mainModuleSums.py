# coding: utf-8
# python3 mainModuleSums.py --event -1 --geom V16 --algo baseline --particle pions --subdet 5

_all_ = [ ]

import os
import sys

parent_dir = os.path.abspath(__file__ + 4 * '/..')
sys.path.insert(0, parent_dir)

import argparse
import algosMS
import processingMS
import resolutionMS
import helperMS
import geometryMS


def parse_arguments():
    parser = argparse.ArgumentParser(description="Interactive Grid Comparison")

    parser.add_argument("--subdet", type=int, default=1, help="1: CEE (has only silicon layers), 2: CEH - only silicon part, 3: CEH - only scint part, 4: CEH, all layers, 5: CEE + CEH")
    parser.add_argument("--event", default='5492', help="Select event number or -1 for all events")
    parser.add_argument("--n", type=int, default=None, help="Process n events (random ordering)")
    parser.add_argument("--algo", default='8towers', help="Select algorithm (baseline, area_overlap, 8towers, 16towers)")
    parser.add_argument("--particle", default='photons', help="Select particle type (photons or pions)")
    parser.add_argument("--geom", default='V16', help="Select the CMSSW geometry (V11 or V16)")
    return parser.parse_args()


def main(subdet, event, particle, algo, n, geom):
    process = processingMS.Processing() 
    algorithms = algosMS.Algorithms()
    resolution = resolutionMS.Resolution()
    helper = helperMS.Helper()
    geometry = geometryMS.Geometry()

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
    elif geom=='V16' and particle == 'pions':
        file_path = (
        '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/'
        'SinglePionPU0V16/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
        )

    # Method that retrieves events and process the data with the geometry
    data, events_to_process = process.get_data_new(event,n,geom,subdet, particle)
    print("DATA", data.columns)
    print("Events_to_process",events_to_process )

    data_gen = process.get_genpart_data(file_path, event, events_to_process, n)

    # Retain only one row per unique event in 'data_gen'
    data_gen = data_gen.drop_duplicates(subset='event', keep='first')
    print("Dataframe columns",data_gen['event'])
    print("DATA GEN", data_gen)

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

    process.ModSumToTowers(initial_kw, data , subdet, event, particle, algo, bin_geojson_filename, hex_geojson_filename, data_gen, geom)

    #geometry.save_bin_geo(towers_bins, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_with_arcs.geojson', f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_only_vertices.geojson')
    #geometry.save_bin_hex(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/hexagons_CMSSW.geojson')
    #geometry.save_scint_mod_geo(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/scint_modules_geo.geojson')

if __name__ == '__main__':
    args = parse_arguments()
    main(args.subdet, args.event, args.particle, args.algo, args.n, args.geom)

    

