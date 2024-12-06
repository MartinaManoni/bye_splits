#python3 produce_OnTier3.py --inputFolder /SinglePion_Pt-0To200-gun/NEW_SinglePion_Pt-0To200-gun_Phase2Fall22DRMiniAOD-noPU_125X_mcRun4_realistic_v2-v1/240615_114023/0000 --outputFolder NEWSinglePion_FlatPt-1To100-gun --n_files 2 --particles pions --queue short



#python3 produce_OnTier3.py --inputFolder /MinBias_TuneCP5_14TeV-pythia8/Min_Bias_try2/241204_145125/0000/ --outputFolder MinBias_TuneCP5_14TeV-pythia8_Phase2Fall22DRMiniAOD-PU200_125X_mcRun4_realistic_v2-v1 --n_files 2 --particles pions --queue short

#root://eos.grif.fr//eos/grif/cms/llr/store/user/mmanoni/MinBias_TuneCP5_14TeV-pythia8/Min_Bias_try2/241204_145125/0000/file_2.root

import os
import glob
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--inputFolder", dest="inputFolder",   type=str, default=None, help="Input folder on the grid")
parser.add_option("--outputFolder",dest="outputFolder",  type=str, default=None, help="Ouutput folder where to store")
parser.add_option("--n_files",     dest="n_files",       type=int, default=5,  help="Number of ntuples to process")
parser.add_option("--particles",   dest="particles",     type=str, default='photons', help="Type of particle to process")
parser.add_option("--queue",       dest="queue",         type=str, default='short', help="long or short queue")
parser.add_option("--no_exec",     dest="no_exec",       action='store_true', default=False, help="stop execution")

#Input Folders ?
#DoublePhoton_FlatPt-1To100-gun/DoublePhoton_FlatPt-1To100-gun_Phase2Fall22DRMiniAOD-PU200_125X_mcRun4_realistic_v2-v1/240613_095947/0000
#SinglePhoton_Pt-2To200-gun/SinglePhoton_Pt-2To200-gun_Phase2Fall22DRMiniAOD-noPU_125X_mcRun4_realistic_v2-v1/240614_081036/0000
#SinglePion_Pt-0To200-gun/SinglePion_Pt-0To200-gun_Phase2Fall22DRMiniAOD-PU200_125X_mcRun4_realistic_v2-v1/240614_140507/0000
#SinglePion_Pt-0To200-gun/SinglePion_Pt-0To200-gun_Phase2Fall22DRMiniAOD-noPU_125X_mcRun4_realistic_v2-v1/240614_135746/0000


#/SinglePion_Pt-0To200-gun/NEW_SinglePion_Pt-0To200-gun_Phase2Fall22DRMiniAOD-noPU_125X_mcRun4_realistic_v2-v1/240615_114023/0000


#Output Folders
#DoublePhoton_FlatPt-1To100-gun
#SinglePhoton_Pt-2To200-gun

(options, args) = parser.parse_args()

infile_base  = os.getcwd()+'/../'
user = infile_base.split('/')[5]
outfile_base = "/data_CMS/cms/"+user+"L1HGCAL/ntupleV16Production/"
infile_base  = "root://eos.grif.fr//eos/grif/cms/llr/store/user/mmanoni/"

###########

folder = outfile_base+options.outputFolder+'/skimmed_ntuples'
queue = options.queue
os.system('mkdir -p ' + folder)

print("Input has" , options.n_files, "files")
for idx in range(options.n_files):
    outRootName = folder + '/Ntuple_' + str(idx+1) + '.root'
    outJobName  = folder + '/job_' + str(idx+1) + '.sh'
    outLogName  = folder + "/log_" + str(idx+1) + ".txt"

    root_file = infile_base+options.inputFolder+"/file_"+str(idx+1)+".root"

    #/SinglePion_Pt-0To200-gun_noPU_125X_mcRun4_realistic_v2-v1_

    #/V16_DoublePhoton_FlatPt-1To100-gun_PU200_125X_mcRun4_realistic_v2-v1_
    #/V16_SinglePhoton_Pt-2To200-gun_noPU_125X_mcRun4_realistic_v2-v1_
    #/SinglePion_Pt-0To200-gun_PU200_125X_mcRun4_realistic_v2-v1_
    #/SinglePion_Pt-0To200-gun_noPU_125X_mcRun4_realistic_v2-v1_
    #SinglePion_Pt-0To200-gun_noPU_125X_mcRun4_realistic_v2-v1_2.root

    #/eos/grif/cms/llr/store/user/mmanoni/SinglePion_Pt-0To200-gun/NEW_SinglePion_Pt-0To200-gun_Phase2Fall22DRMiniAOD-noPU_125X_mcRun4_realistic_v2-v1/240615_114023
    cmsRun = "python3 produce_mod_MinBias.py --nevents=-1 --particles "+options.particles+" --inputFile="+root_file+" --outputFile="+outRootName+" >& "+outLogName

    #produce_mod_updated.py
    skimjob = open (outJobName, 'w')
    skimjob.write ('#!/bin/bash\n')
    skimjob.write ('export X509_USER_PROXY=~/.t3/proxy.cert\n')
    skimjob.write ('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    skimjob.write ('cd %s\n' % os.getcwd())
    #skimjob.write ('export SCRAM_ARCH=slc6_amd64_gcc472\n')
    skimjob.write ('eval `scram r -sh`\n')
    skimjob.write ('cd %s\n'%os.getcwd())
    skimjob.write (cmsRun+'\n')
    skimjob.close ()

    os.system ('chmod u+rwx ' + outJobName)
    command = ('/home/llr/cms/'+user+'/t3submit -'+queue+' \'' + outJobName +"\'")
    print(command)
    if not options.no_exec: os.system (command)
