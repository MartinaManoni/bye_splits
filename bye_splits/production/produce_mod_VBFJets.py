# coding: utf-8
# python3 bye_splits/production/produce_mod.py --particles photons --nevents 3000
_all_ = []

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

from bye_splits.utils import params

import ROOT
import yaml

ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<float> calcDeltaR(ROOT::VecOps::RVec<float> geta, ROOT::VecOps::RVec<float> gphi,
                                     ROOT::VecOps::RVec<float> cleta, ROOT::VecOps::RVec<float> clphi) {
    if (geta.size() == 0) { // empty event (filtered before)
        return ROOT::VecOps::RVec<float>();
    }
    assert(geta.size() == 1); // consistency check
    unsigned ncl = cleta.size();
    ROOT::VecOps::RVec<float> deltaR(ncl);
    float deta, dphi;

    for (unsigned j = 0; j < ncl; ++j) {
        deta = fabs(cleta[j] - geta[0]);
        dphi = fabs(clphi[j] - gphi[0]);
        if (dphi > M_PI) dphi -= (2 * M_PI);
        deltaR[j] = sqrtf(dphi * dphi + deta * deta);
    }

    return deltaR;
}
""")

ROOT.gInterpreter.Declare("""
    ROOT::RDF::RResultPtr<ULong64_t> addProgressBar(ROOT::RDF::RNode df) {
        auto c = df.Count();
        c.OnPartialResult(/*every=*/100, [](ULong64_t e) { std::cout << e << std::endl; });
        return c;
    }
""")

ROOT.gInterpreter.Declare("""
vector<int> convertInt(const ROOT::VecOps::RVec<int> &v) {
    return vector<int>(v.begin(), v.end());
}
""")
ROOT.gInterpreter.Declare("""
vector<unsigned> convertUint(const ROOT::VecOps::RVec<unsigned> &v) {
    return vector<unsigned>(v.begin(), v.end());
}
""")
ROOT.gInterpreter.Declare("""
vector<float> convertFloat(const ROOT::VecOps::RVec<float> &v) {
    return vector<float>(v.begin(), v.end());
}
""")
ROOT.gInterpreter.Declare("""
vector<vector<float>> convertFloat2D(const ROOT::VecOps::RVec<vector<float>> &v) {
    vector<vector<float>> vec(v.size());
    for (unsigned i = 0; i < v.size(); ++i) {
        vec[i] = vector<float>(v[i].begin(), v[i].end());
    }
    return vec;
}
""")

ROOT.gInterpreter.Declare("""
#include <set>
#include <cmath>
#include <vector>
#include <iostream> // for print statements

std::tuple<ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<int>> GetVBFEventFlags(
    const ROOT::VecOps::RVec<int>& particle_pdgid,
    const ROOT::VecOps::RVec<int>& particle_status,
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& particle_daughters,
    const ROOT::VecOps::RVec<float>& particle_eta,
    const ROOT::VecOps::RVec<float>& particle_phi,
    const ROOT::VecOps::RVec<float>& jet_eta,
    const ROOT::VecOps::RVec<float>& jet_phi
) {
    // Threshold for delta R
    const float deltaR_threshold = 0.1;

    // To store indices of processed daughters
    std::set<int> processed_daughters;

    // Initialize flags for jets (all set to false initially)
    ROOT::VecOps::RVec<bool> jet_flags(jet_eta.size(), false);

    // Initialize event flag to good (0 by default)
    int event_flag = 0; // 0 means good event, 1 means bad event

    // Iterate over all particles
    for (size_t i = 0; i < particle_pdgid.size(); ++i) {
        int particle_pdgid_current = particle_pdgid[i];
        int particle_status_current = particle_status[i];

        // Check if the particle is a gluon (PDG ID 21) with status 21
        if (particle_pdgid_current == 21 && particle_status_current == 21) {
            // If a gluon mother with status 21 is found, flag the event as bad
            event_flag = 1; // Mark event as bad
            std::cout << "Bad event: Found gluon (PDG ID: 21) with status 21" << std::endl;
            break; // No need to process further if the event is bad
        }
    }

    // Only proceed to quark mother processing if the event is not bad
    if (event_flag == 0) {
        // Iterate over all particles again to check for quark mothers with status 21
        for (size_t i = 0; i < particle_pdgid.size(); ++i) {
            int particle_pdgid_current = particle_pdgid[i];
            int particle_status_current = particle_status[i];

            // Check if the particle is a quark (PDG ID 1-6 or -1 to -6) with status 21
            if (std::abs(particle_pdgid_current) >= 1 && std::abs(particle_pdgid_current) <= 6 && particle_status_current == 21) {
                // Print out the particle information
                //std::cout << "Found quark mother (PDG ID: " << particle_pdgid_current 
                        //  << ") with status 21. Checking daughters..." << std::endl;

                // Retrieve the daughters of the current particle
                const ROOT::VecOps::RVec<int>& daughters_indices = particle_daughters[i];

                // Process the quark daughters
                for (auto daughter_index : daughters_indices) {
                    if (daughter_index >= 0 && daughter_index < static_cast<int>(particle_pdgid.size())) {
                        // Skip if the daughter has already been processed
                        if (processed_daughters.find(daughter_index) != processed_daughters.end()) {
                            continue;
                        }

                        int daughter_pdgid = particle_pdgid[daughter_index];
                        if (std::abs(daughter_pdgid) >= 1 && std::abs(daughter_pdgid) <= 6) { // Check if it's a quark
                            processed_daughters.insert(daughter_index);

                            float eta = particle_eta[daughter_index];
                            float phi = particle_phi[daughter_index];

                            // Print out daughter particle details
                            //std::cout << "Processing quark daughter (PDG ID: " << daughter_pdgid
                                    //  << ") with eta: " << eta << " and phi: " << phi << std::endl;

                            // Check matching with gen jets
                            for (size_t j = 0; j < jet_eta.size(); ++j) {
                                float delta_eta = eta - jet_eta[j];
                                float delta_phi = std::abs(phi - jet_phi[j]);
                                if (delta_phi > M_PI) delta_phi = 2 * M_PI - delta_phi; // Wrap-around for phi
                                float deltaR = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);

                                if (deltaR < deltaR_threshold) {
                                    jet_flags[j] = true; // Set flag to true for matched jet
                                    //std::cout << "Match found for quark daughter (PDG ID: " << daughter_pdgid
                                              //<< ") with gen jet | ΔR: " << deltaR << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // If the event is bad, mark all jet flags as false (no need to check jets)
    if (event_flag == 1) {
        for (size_t j = 0; j < jet_flags.size(); ++j) {
            jet_flags[j] = false; // Mark the jets as unmatched for a bad event
        }
    }

    // Print event flag status
    //std::cout << "Event flag: " << event_flag << " (0: good, 1: bad)" << std::endl;

    // Return separate event flag and jet flags
    ROOT::VecOps::RVec<int> output_event_flag = {event_flag};
    ROOT::VecOps::RVec<int> output_jet_flags;
    for (size_t j = 0; j < jet_flags.size(); ++j) {
        output_jet_flags.push_back(jet_flags[j] ? 1 : 0); // Convert jet flags to 1/0
    }

    // Return both event flag and jet flags
    return std::make_tuple(output_event_flag, output_jet_flags);
}
""")




# The skim function takes several parameters,
# including the tree name (tn), input file (inf), output file (outf),
# particle type (particle), number of events to process (nevents),
# and a configuration dictionary (cfg)

def skim(tn, inf, outf, particle, nevents, cfg):
    if cfg["selection"]["disconnectedTriggerLayers"]:
        discLayers = cfg["selection"]["disconnectedTriggerLayers"]
    if cfg["selection"]["reachedEE"]:
        reachedEE = cfg["selection"]["reachedEE"]
    if cfg["selection"]["deltarThreshold"]:
        deltarThreshold = cfg["selection"]["deltarThreshold"]
    if cfg["selection"]["mipThreshold"]:
        mipThreshold = cfg["selection"]["mipThreshold"]

    if nevents == -1:  # RDataFrame.Range() does not work with multithreading
        ROOT.EnableImplicitMT()
        print("Multithreaded...")

    dataframe = ROOT.RDataFrame(tn, inf)

    # gen related variables
    gen_intv = ["gen_pdgid", "gen_charge", "gen_status", "gen_daughters"]
    gen_floatv = ["gen_eta", "gen_phi", "gen_pt", "gen_energy"]
    gen_v = gen_intv + gen_floatv
    print("gen_v", gen_v)

    #dd=dataframe.Define("tmp_good_gens") #no conditions applied for gen particles
    condgen = "gen_pt >= 0"
    dd = dataframe.Define("tmp_good_gens", condgen)
    for v in gen_v:
        print("v", v)
        dd = dd.Define("tmp_good_" + v, v + "[tmp_good_gens]")
    
    # jets related variables
    #genjet_intv = ["genjet_n"]
    genjet_floatv = ["genjet_pt", "genjet_energy", "genjet_eta", "genjet_phi"]
    genjet_v= genjet_floatv

    condgenjets = "genjet_pt >= 0"
    dd0 = dd.Define("tmp_good_genjets", condgenjets)
    for v in genjet_v:
        dd0 = dd0.Define("tmp_good_" + v, v + "[tmp_good_genjets]")

    print(dd0.GetColumnNames())

    # trigger cells-related variables
    tc_uintv = ["tc_multicluster_id"]
    tc_intv = ["tc_layer", "tc_cellu", "tc_cellv", "tc_waferu", "tc_waferv", "tc_subdet"]
    tc_floatv = ["tc_energy", "tc_mipPt", "tc_pt", "tc_x", "tc_y", "tc_z", "tc_phi", "tc_eta"]
    tc_v = tc_uintv + tc_intv + tc_floatv


    # selection on trigger cells (within each event)
    condtc = "tc_mipPt > " + str(mipThreshold)
    dd1 = dd0.Define("tmp_good_tcs", condtc)
    for v in tc_v:
        dd1 = dd1.Define("tmp_good_" + v, v + "[tmp_good_tcs]")
    
    # trigger sums-related variables
    ts_intv = ["ts_layer", "ts_waferu", "ts_waferv", "ts_subdet"]
    ts_floatv = ["ts_energy", "ts_mipPt", "ts_pt", "ts_x", "ts_y", "ts_z", "ts_phi", "ts_eta"]
    ts_v = ts_intv + ts_floatv

    # selection on trigger sums (within each event)
    condts = "ts_mipPt > 0"
    dd1 = dd1.Define("tmp_good_ts", condts)
    for v in ts_v:
        dd1 = dd1.Define("tmp_good_" + v, v + "[tmp_good_ts]")  

    # cluster-related variables
    cl_uintv = ["cl3d_id"]
    cl_floatv = ["cl3d_energy", "cl3d_pt", "cl3d_eta", "cl3d_phi"]
    cl_v = cl_uintv + cl_floatv

    # selection on clusters (within each event)
    condcl = "cl3d_eta > 0"
    dd1 = dd1.Define("tmp_good_cl", condcl)
    for v in cl_v:
        dd1 = dd1.Define("tmp_good_" + v, v + "[tmp_good_cl]")

    # remove events with zero clusters
    dd2 = dd1.Filter("tmp_good_cl3d_id.size()!=0")


    dd2 = dd2.Define("event_flag_jet_flags", 
    "GetVBFEventFlags(tmp_good_gen_pdgid, tmp_good_gen_status, tmp_good_gen_daughters,tmp_good_gen_eta, tmp_good_gen_phi, tmp_good_genjet_eta, tmp_good_genjet_phi)")

    #Split the tuple into separate columns for event flag and jet flags
    matchvars = ["event_flag", "jet_flag"]
    dd2 = dd2.Define(matchvars[0], "convertInt(std::get<0>(event_flag_jet_flags))").Define(matchvars[1], "convertInt(std::get<1>(event_flag_jet_flags))")

    print(dd2.GetColumnNames())

    # convert root vector types to vector equivalents (uproot friendly)
    intv = tc_intv + ts_intv + gen_intv
    for var in intv:
        if "gen_daughters" in var:
            # Flatten the gen_daughters (vector of vectors) into a 1D vector
            dd2 = dd2.Define("good_" + var, 
                            "std::vector<int> flat_daughters; \
                            for (size_t i = 0; i < tmp_good_" + var + ".size(); ++i) { \
                                for (size_t j = 0; j < tmp_good_" + var + "[i].size(); ++j) { \
                                    flat_daughters.push_back(tmp_good_" + var + "[i][j]); \
                                } \
                            } \
                            return flat_daughters;")
        else:
            dd2 = dd2.Define("good_" + var, "convertInt(tmp_good_" + var + ")")

    uintv = cl_uintv + tc_uintv
    for var in uintv:
        dd2 = dd2.Define("good_" + var, "convertUint(tmp_good_" + var + ")")

    floatv = tc_floatv + cl_floatv + ts_floatv + gen_floatv+ genjet_floatv
    for var in floatv:
        dd2 = dd2.Define("good_" + var, "convertFloat(tmp_good_" + var + ")")

    # define stored variables (and rename some)
    allvars = tc_v + cl_v + ts_v + gen_v + genjet_v
    good_allvars = ["event"] + matchvars
    for v in allvars:
        print("all vars", v)
        good_allvars.append("good_" + v)
    for v in good_allvars:
        print("GOOD all vars", v)
    # store skimmed file
    if nevents == -1:
        dd2.Snapshot(tn, outf, good_allvars)
    elif nevents > 0:
        dd2.Range(0, nevents).Snapshot(tn, outf, good_allvars)
    else:
        dd2.Snapshot(tn, outf, good_allvars)

    print('OutputFolder:', outf)
    #dd2.Count().OnPartialResult(10, "[&log](auto c) { l << c << \" events processed\n\";}")

    # display event processing progress
    count = ROOT.addProgressBar(ROOT.RDF.AsRNode(dd2))
    count.GetValue()


if __name__ == "__main__":
    ''' produce.py <nevents> <inputFile> <outputFile> '''
    import argparse
    parser = argparse.ArgumentParser(description='Skim ntuples')
    parser.add_argument('--particles', type=str, choices=('photons', 'electrons', 'pions'),
                        required=False, default='photons', help='particles to skim')
    parser.add_argument('--nevents', type=int, default=100,
                        required=False, help='number of events to skim')
    parser.add_argument('--inputFile', type=str, default=None,
                        required=False, help='root file in input')
    parser.add_argument('--outputFile', type=str, default=None,
                        required=False, help='skimmed root file in output')
    FLAGS = parser.parse_args()

    with open(params.CfgPath, 'r') as afile:
        cfg = yaml.safe_load(afile)

    dir_tree = os.path.join(cfg["io"]["production"][FLAGS.particles]["dir"],
                            cfg["io"]["production"][FLAGS.particles]["tree"])
    
    if FLAGS.inputFile and FLAGS.outputFile:
        infile, outfile = FLAGS.inputFile, FLAGS.outputFile
    else:
        infile = cfg["io"]["production"][FLAGS.particles]["infile"]
        outfile = cfg["io"]["production"][FLAGS.particles]["outfile"]
    skim(dir_tree,"/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/production/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8_Phase2Spring23DIGIRECOMiniAOD.root" ,"/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/production/output_VBF.root" , FLAGS.particles, FLAGS.nevents, cfg)
