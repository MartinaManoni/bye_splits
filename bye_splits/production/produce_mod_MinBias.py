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

    # trigger cells-related variables
    tc_uintv = ["tc_multicluster_id"]
    tc_intv = ["tc_layer", "tc_cellu", "tc_cellv", "tc_waferu", "tc_waferv", "tc_subdet"]
    tc_floatv = ["tc_energy", "tc_mipPt", "tc_pt", "tc_x", "tc_y", "tc_z", "tc_phi", "tc_eta"]
    tc_v = tc_uintv + tc_intv + tc_floatv


    # selection on trigger cells (within each event)
    condtc = "tc_mipPt > " + str(mipThreshold)
    dd1 = dataframe.Define("tmp_good_tcs", condtc)
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

    # convert root vector types to vector equivalents (uproot friendly)
    intv = tc_intv + ts_intv
    for var in intv:
        dd2 = dd2.Define("good_" + var, "convertInt(tmp_good_" + var + ")")

    uintv = cl_uintv + tc_uintv
    for var in uintv:
        dd2 = dd2.Define("good_" + var, "convertUint(tmp_good_" + var + ")")

    floatv = tc_floatv + cl_floatv + ts_floatv
    for var in floatv:
        dd2 = dd2.Define("good_" + var, "convertFloat(tmp_good_" + var + ")")

    # define stored variables (and rename some)
    allvars = tc_v + cl_v + ts_v
    good_allvars = ["event"]
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
    skim(dir_tree, infile, outfile, FLAGS.particles, FLAGS.nevents, cfg)
