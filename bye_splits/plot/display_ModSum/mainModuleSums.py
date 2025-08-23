# coding: utf-8
# python3 mainModuleSums.py --event -1 --geom V16 --algo baseline --particle pions --subdet 5
# python3 mainModuleSums.py --event -1 --n 10 --geom V16 --algo 8towers --particle neutrinos --subdet 5 --inputfile root --STCs
# python3 mainModuleSums.py --event -1 --n 4000 --geom V16 --algo 4towers --particle pions --subdet 5 --inputfile root --STCs

# python3 mainModuleSums.py --event -1 --n 4000 --geom V16 --algo 4towers --particle photons --subdet 5 --inputfile root --no-STCs
_all_ = [ ]

import os
import sys
import re

parent_dir = os.path.abspath(__file__ + 4 * '/..')
sys.path.insert(0, parent_dir)

import argparse
import processingMS
import resolutionMS
import helperMS
import geometryMS
import json
import warnings
import pandas as pd
warnings.filterwarnings("ignore", category=UserWarning, message=".*subnormal.*")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Interactive Grid Comparison")

    parser.add_argument("--subdet", type=int, default=1, help="1: CEE (has only silicon layers), 2: CEH - only silicon part, 3: CEH - only scint part, 4: CEH, all layers, 5: CEE + CEH")
    parser.add_argument("--event", default='5492', help="Select event number or -1 for all events")
    parser.add_argument("--n", type=int, default=None, help="Process n events (random ordering)")
    parser.add_argument("--algo", default='8towers', help="Select algorithm (baseline, area_overlap, 4towers, 8towers, 16towers)")
    parser.add_argument("--particle", default='photons', help="Select particle type (photons, pions or neutrinos)")
    parser.add_argument("--geom", default='V16', help="Select the CMSSW geometry (V11 or V16)")
    parser.add_argument("--inputfile", default='root', help="Select input file type (root or hdf5)")
    parser.add_argument("--inputdir", default=None, help="Path to directory with multiple ROOT files")
    parser.add_argument("--outputdir", default="./outputs", help="Directory to store output files")
    parser.add_argument('--range', nargs=2, type=int, metavar=('START', 'END'),
                    help='Range of ntuples to process (e.g., --range 1 20)')

    # Boolean flag to toggle STCs option
    parser.add_argument("--PU200", dest='PU200', action='store_true', help="Use PU200 sample for VBF")
    parser.add_argument("--STCs", dest='STCs', action='store_true', help="Enable STCs")
    parser.add_argument("--no-STCs", dest='STCs', action='store_false', help="Disable STCs")
    return parser.parse_args()


def main(subdet, event, particle, algo, n, geom, inputfile, STCs, PU200, root_file=None, output_file=None):
    process = processingMS.Processing()
    geometry = geometryMS.GeometryData()
    print("Entering Main Function ....")

    save_tower_bins = False

    #file_path = f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/DoublePhotonsPU0_3k_V11/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
    if root_file is None:
        if geom=='V11':
            print("using V11 geom")
            hdf5_file = (
            '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/'
            'DoublePhotonsPU0_hadd_123_energy/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
            )
        elif geom=='V16' and particle == 'photons':
            hdf5_file = (
            '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/'
            'SinglePhotonPU0V16/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
            )
            root_file =('/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples/SinglePhotonPU0V16.root')

        elif geom=='V16' and particle == 'pions':
            hdf5_file = (
            '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/'
            'SinglePionPU0V16/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5'
            )
            root_file =('/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples/SinglePionPU0V16.root')
            if STCs:
                print("Enabling STCs for Pions samples...")
                root_file =('/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples_STCS/SinglePionPU0V16_STCs.root')

        elif geom=='V16' and particle == 'neutrinos':
            root_file =('/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples/MinBiasPU200_Fall22.root')
            if STCs:
                print("Enabling STCs for MinBias samples...")
                root_file =('/data_CMS/cms/manoniL1HGCAL/ntupleV16Production/MinBias_STCs_Final/skimmed_ntuples/Ntuple_1.root')

        elif geom == 'V16' and particle == 'jets':
            if root_file is not None:
                print(f"Using provided ROOT file: {root_file}")
            elif STCs and PU200:
                print("Enabling STCs for Jets samples PU200...")
                root_file = '/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples_STCS/VBFHToInvisible_Spring23_PU200_100ntuples.root'
            else:
                print("Enabling STCs for Jets samples PU0...")
                root_file = '/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples/VBFHtoInvPU0_STCs_Ntuple12.root'
    else:
        print(f"Using provided ROOT file: {root_file}")

    if inputfile == "hdf5":
        print("processing hdf5 file....")
        data, events_to_process = process.get_data_new(event, n, geom, subdet, particle)
        df_gen, events = process.get_genpart_data(hdf5_file, event, events_to_process, n)
        df_specific = data

    else:
        print("processing root file....")
        # Skip get_gen_particles for neutrinos
        if particle == 'neutrinos':
            print("Skipping get_gen_particles for neutrinos")
            df_gen, events = None, None  # Set these to None since gen variables are unavailable
        else:
            df_gen, events = process.get_gen_particles(root_file, particle, n, event)
            print("df_gen col", df_gen.columns)
            print("events", events)

        # Create the specific data frame for the particle
        if events is None or len(events) == 0:
            selected_events = None
        else:
            selected_events = events


        #print("selected_events", len(selected_events))

        if STCs:
            if subdet == 1:
                # Process df_specific only for subdet 1
                print('processing subdet 1 - STCs')
                df_specific = process.read_root_and_create_dataframe(
                    root_file, subdet, selected_events=selected_events
                )
            elif subdet in [2, 3]: #per ora implemento solo singoli subdet no somma
                # Process df_STCs for subdet 2, 3
                #read_root_and_create_dataframe_STCS non contine process_V16 perche non necessario, solo scelgi il subet adatto
                print(f'processing {subdet} - STCs')
                df_STCs = process.read_root_and_create_dataframe_STCS(
                    root_file, subdet,selected_events=selected_events
                )
                df_specific = None
            elif subdet in [5]:
                print(f'processing {subdet} - STCs')
                df_specific = process.read_root_and_create_dataframe(
                    root_file, 1, selected_events=selected_events
                )
                df_STCs_2 = process.read_root_and_create_dataframe_STCS(
                    root_file, 2,selected_events=selected_events
                )

                df_STCs_3 = process.read_root_and_create_dataframe_STCS(
                    root_file, 3,selected_events=selected_events
                )

                # Merge df_STCs_2 and df_STCs_3
                df_STCs = pd.concat([df_STCs_2, df_STCs_3], ignore_index=True)
            else:
                raise ValueError(f"Unsupported subdet value: {subdet}")

        else:
            df_specific = process.read_root_and_create_dataframe(
            root_file, subdet, selected_events= selected_events)
            df_STCs = None

        print("data col", df_specific.columns)

    #helper.read_hdf5_structure(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/photons_manoni/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5')
    #helper.read_all_block0_values(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/photons_manoni/fill_gencl_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5')

    bin_geojson_filename = '/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_with_arcs.geojson' #bins_with_arcs
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
        print("creating and saving tower bins...")
        process.create_and_save_tower_bins(initial_kw, df_specific, geom) #create and save tower bins

    process.ModSumToTowers(initial_kw, df_specific, df_STCs, subdet, event, particle, algo, bin_geojson_filename, hex_geojson_filename, df_gen, geom, STCs, output_file)

    #geometry.save_bin_geo(towers_bins, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_with_arcs.geojson', f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_only_vertices.geojson')
    #geometry.save_bin_hex(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/hexagons_CMSSW.geojson')
    #geometry.save_scint_mod_geo(f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/scint_modules_geo.geojson')

if __name__ == '__main__':
    args = parse_arguments()

    def natural_sort_key(s):
    # Extract numbers and strings to sort naturally
        return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

    if args.inputdir is not None:
        import glob
        root_files = glob.glob(os.path.join(args.inputdir, "*.root"))
        root_files = sorted(root_files, key=natural_sort_key)
        print(f"Found {len(root_files)} ROOT files in inputdir.")

        # Apply --range if provided
        if args.range:
            start_idx, end_idx = args.range
            root_files = root_files[start_idx - 1:end_idx]  # Convert to 0-based index
            print(f"Processing files from Ntuple_{start_idx} to Ntuple_{end_idx}")
        else:
            root_files = root_files[:100]

        for idx, rf in enumerate(root_files):
            print(f"\nProcessing file {idx+1}/{len(root_files)}: {rf}")

            base_name = os.path.splitext(os.path.basename(rf))[0] #Ntuple_1
            output_filename = f"{base_name}_output.txt" 
            print(f"Will save output to: {output_filename}")

            # Copy args to modify safely, disable inputdir to avoid recursion
            args_copy = argparse.Namespace(**vars(args))
            args_copy.inputdir = None
            args_copy.inputfile = 'root'

            # Pass current root file to main()
            main(args_copy.subdet, args_copy.event, args_copy.particle, args_copy.algo,
                 args_copy.n, args_copy.geom, args_copy.inputfile, args_copy.STCs, args_copy.PU200,
                 root_file=rf,output_file=output_filename)
    else:
        main(args.subdet, args.event, args.particle, args.algo, args.n, args.geom, args.inputfile, args.STCs, args.PU200)

    

