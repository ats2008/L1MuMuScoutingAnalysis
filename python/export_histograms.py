import uproot as urt
import json,os
import numpy as np
import awkward as ak
import vector,argparse
vector.register_awkward()

from util import *

base='workarea/datasets'
file_list={
     'singleMu'   : base+'/15_0_1_pre2/singleMu_perfNano.root',
    'DY1Jet_Mll_1To10' : base+'/DY1Jet_Mll_1To10.root',
    'DY_Mll_10To50'    : base+'/DY_Mll_10To50.root',
    'MinBias'          : base+'/MinBias.root',
    'DPhoton_M1'       : base+'/DARK_PHOTON_M1.root',
    'DPhoton_M2'       : base+'/DARK_PHOTON_M2.root',
    'DPhoton_M4'       : base+'/DARK_PHOTON_M4.root',
    'DPhoton_M8'       : base+'/DARK_PHOTON_M8.root',
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inputFile", help="inputFile")
    parser.add_argument("-t","--tag", help=" tag")
    parser.add_argument("-d","--destination", help="destination folder",default='./')
    parser.add_argument("-r","--doc", help="documentation string",default='MuMuSel')
    parser.add_argument("-e","--eMax", help="max Events",default=-1,type=int)
    args = parser.parse_args()
    tag=args.tag
    prefix=args.destination
    file_list={
        tag : args.inputFile,
    }
    treeStore={}
    fileToLoad=list(file_list.keys())
    for tag in fileToLoad:
        fl=file_list[tag]
        print(f"Loading {tag} : {fl}")
        file=urt.open(fl)
        tree=file['Events']
        treeStore[tag]=tree
    #     break
        print(f"   the event tree has {treeStore[tag].num_entries} events")
    
    allBranches=[]
    _set=[]
    for ky in tree.keys():
    #     print(ky)
        skip=True
        if 'TkMuon' in ky:
            skip=False
        if 'GenCands' in ky:
            skip=False
        if 'L1PFCands' in ky:
            skip=False
        if 'L1PuppiCands' in ky:
            skip=False
        if 'L1Vertex' in ky:
            skip=False
        if skip : continue
        allBranches.append(ky)
        if ky.split("_")[0] not in _set:
            _set.append( ky.split("_")[0] )   
            if not _set[-1].startswith('n'):
                print("Adding banches : ",_set[-1])
    dataStore={}
    for tag in treeStore:
        print("Loading arrays for ",tag)
        if args.eMax>0:
            dataStore[tag]=treeStore[tag].arrays(allBranches,entry_stop=args.eMax)
        else:
            dataStore[tag]=treeStore[tag].arrays(allBranches)
        print(f"   loaded {len(dataStore[tag])} events")
    
    data=dataStore[tag]
    
    config_dict={}
    config_dict["METADATA"]={}
    config_dict["METADATA"]["tag"]=tag
    config_dict["METADATA"]["filename"]=file_list[tag]
    config_dict["METADATA"]["NEVENTS"]=len(dataStore[tag])
    config_dict["METADATA"]["documentation"]=args.doc
    config_dict["NEVENTS"]=len(dataStore[tag])
    
    config_dict["SELECTION_COUNTS"]={}
    
    config_dict["HISTSTORE"]={}
    
    histStore=config_dict["HISTSTORE"]
    histStore["nTkMuons"]=getHistograms(data.nL1gmtTkMuon,np.arange(0.0,20,1)-0.5,"nTkMuons","Number of tk muon candidates prior to any selection")
    histStore["nPF"]     =getHistograms(data.nL1PFCands,np.arange(150.0,750,5),"nPF","Number of pf candidates prior to any selection")
    histStore["nPuppi"]  =getHistograms(data.nL1PuppiCands,np.arange(0.0,100,5),"nPuppi","Number of puppi candidates prior to any selection")
    
    dataStore=addGenMatchedMuonInformation(dataStore,pTMinMuon=2.0,etaMax=2.4)
    #dataStore=addGenMatchedMuonInformation(dataStore,pTMinMuon=4.0,etaMax=1.9)
    
    ARRAY=ak.flatten(dataStore[tag]['matchedDimuons'].mass)
    BINS=np.arange(0.0,20.0,0.1)
    NAME="genMatchedDimuMass"
    DOC="Gen matched tk di-muons"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
    
    
    ARRAY=dataStore[tag]['matchedDimuons'].pt
    BINS=np.arange(0.0,40.0,0.5)
    NAME="genMatchedDimuPt"
    DOC="Gen matched tk di-muons pt"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
    
    dataStore=addSelectedDimuons(dataStore,pTMinMuon=2.0,etaMax=2.4,d0Max=0.05)
    #dataStore=addSelectedDimuons(dataStore,pTMinMuon=2.0,etaMax=2.4)
    #dataStore=addSelectedDimuons(dataStore,pTMinMuon=4.0,etaMax=1.9)
    
    ARRAY=ak.num(dataStore[tag]['selected_dimu'])
    BINS=np.arange(0.0,7.0,1)-0.5
    NAME="dimuonMultiplicityPostSelection"
    DOC="number of dimuon candidate per event, after selection"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
    
    ARRAY=dataStore[tag]['selected_dimu'].mass
    BINS=np.arange(0.0,60.0,0.25)
    NAME="dimuonMass"
    DOC="mass of dimuon candidates, after selection"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
 
    ARRAY=dataStore[tag]['selected_dimu'].mass
    BINS=np.arange(0.0,60.0,0.05)
    NAME="dimuonMassFineBin50Mev"
    DOC="mass of dimuon candidates, after selection, in 50 MeV bins"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

    ARRAY=dataStore[tag]['selected_dimu'].pt
    BINS=np.arange(0.0,60.0,0.25)
    NAME="dimuonPt"
    DOC="pt of dimuon candidates, after selection"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
    
    ARRAY=dataStore[tag]['selected_dimu_mu1'].deltaR(dataStore[tag]['selected_dimu_mu2'])
    BINS=np.arange(0.0,4.0,0.1)
    NAME="dimuonsDrMuMu"
    DOC="deltaR between nuons of dimuon candidates, after selection"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
    
    ARRAY=dataStore[tag]['selected_dimu_mu1'].d0
    BINS=np.arange(0.0,4.0,0.02)
    NAME="dimuonsMu1_d0"
    DOC="d0 / mu1 of  dimuon candidates, after selection"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

    ARRAY=dataStore[tag]['selected_dimu_mu2'].d0
    BINS=np.arange(0.0,4.0,0.02)
    NAME="dimuonsMu2_d0"
    DOC="d0 / mu2 of  dimuon candidates, after selection"
    histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

    ofname=f"{prefix}/out_data_{tag}.json"   
    os.system(f"mkdir  -p {prefix}")
    print(f"exporting to {ofname}")

    with open(ofname,"w") as f:
        json.dump(config_dict,f,indent=4)

if __name__=='__main__':
    main()

