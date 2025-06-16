import uproot as urt
import pickle as pkl
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

branchesToStore=['selected_dimu','selected_dimu_mu1','selected_dimu_mu2']

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
    #dataStore=addSelectedDimuons(dataStore,pTMinMuon=2.0,etaMax=2.4,addTkIsoMask=False,puppiIsoThr=None)
    dataStore=addSelectedDimuons(dataStore,pTMinMuon=2.0,etaMax=2.4,addTkIsoMask=False,puppiIsoThr=0.2001)
    #dataStore=addGenMatchedMuonInformation(dataStore,pTMinMuon=4.0,etaMax=1.9)
    #dataStore=addSelectedDimuons(dataStore,pTMinMuon=4.0,etaMax=1.9)
    
    dataStore=addPFIsolations(dataStore)
    dataStore=addPuppiIsolation(dataStore)

    cat_config={}
    cat_config['highptMuMuCentral']={}
    cat_config['highptMuMuCentral']['min_pt_mu1']=4.0
    cat_config['highptMuMuCentral']['min_pt_mu2']=4.0
    cat_config['highptMuMuCentral']['max_eta_mu1']=1.9
    cat_config['highptMuMuCentral']['max_eta_mu2']=1.9
    
    cat_config['highptMuMuForward']={}
    cat_config['highptMuMuForward']['min_pt_mu1']=4.0
    cat_config['highptMuMuForward']['min_pt_mu2']=4.0
    cat_config['highptMuMuForward']['min_eta_mu1']=1.9
    cat_config['highptMuMuForward']['min_eta_mu2']=1.9
    cat_config['highptMuMuForward']['max_eta_mu1']=2.5
    cat_config['highptMuMuForward']['max_eta_mu2']=2.5
    
    cat_config['MuMuCentral']={}
    cat_config['MuMuCentral']['max_eta_mu1']=1.9
    cat_config['MuMuCentral']['max_eta_mu2']=1.9
    
 
    cat_config['MuMuForeward']={}
    cat_config['MuMuForeward']['min_eta_mu1']=1.9
    cat_config['MuMuForeward']['min_eta_mu2']=1.9
    cat_config['MuMuForeward']['max_eta_mu1']=2.5
    cat_config['MuMuForeward']['max_eta_mu2']=2.5
    
 
    
    categorizedData={}
    categorizedData['inclusive']=dataStore[tag]
    current=categorizedData['inclusive']
    for ky in cat_config:
        current,rest=getCategories( current ,cat_config[ky] )
        categorizedData[ky]=current
        print("CAT : ",ky)
        print("   > Min pt " ,ak.min( current['selected_dimu_mu1'].pt ) ," / ", ak.min( current['selected_dimu_mu2'].pt )   )
        print("   > Max pt " ,ak.max( current['selected_dimu_mu1'].pt ) ," / ", ak.max( current['selected_dimu_mu2'].pt )   )
        print("   > Min eta ",ak.min( np.abs(current['selected_dimu_mu1'].eta) ) ," / ", ak.min( np.abs(current['selected_dimu_mu2'].eta ) )  )
        print("   > Max eta ",ak.max( np.abs(current['selected_dimu_mu1'].eta) ) ," / ", ak.max( np.abs(current['selected_dimu_mu2'].eta ) )  )
        current=rest
    categorizedData['rest']=rest
    
    for ky in list(categorizedData.keys()):
        print(f"Processing {ky} with {len(categorizedData[ky])} events")
        data=categorizedData[ky]
        data=cleanEmptyDimuons(data,cleanEvents=True)
        print(f"    After cleaning  {len(data)} events remains")

        ARRAY=ak.flatten(data['matchedDimuons'].mass)
        BINS=np.arange(0.0,20.0,0.1)
        NAME="genMatchedDimuMass"
        DOC="Gen matched tk di-muons"
        histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
        
        ARRAY=data['matchedDimuons'].pt
        BINS=np.arange(0.0,40.0,0.5)
        NAME="genMatchedDimuPt"
        DOC="Gen matched tk di-muons pt"
        histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
        
        
        ARRAY=ak.num(data['selected_dimu'])
        BINS=np.arange(0.0,7.0,1)-0.5
        NAME="dimuonMultiplicityPostSelection"
        DOC="number of dimuon candidate per event, after selection"
        histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
        
        ARRAY=data['selected_dimu'].mass
        BINS=np.arange(0.0,60.0,0.25)
        NAME="dimuonMass"
        DOC="mass of dimuon candidates, after selection"
        histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
     
        ARRAY=data['selected_dimu'].mass
        BINS=np.arange(0.0,60.0,0.05)
        NAME="dimuonMassFineBin50Mev"
        DOC="mass of dimuon candidates, after selection, in 50 MeV bins"
        histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
    
        ARRAY=data['selected_dimu'].pt
        BINS=np.arange(0.0,60.0,0.25)
        NAME="dimuonPt"
        DOC="pt of dimuon candidates, after selection"
        histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
        
        ARRAY=data['selected_dimu_mu1'].deltaR(data['selected_dimu_mu2'])
        BINS=np.arange(0.0,4.0,0.1)
        NAME="dimuonsDrMuMu"
        DOC="deltaR between nuons of dimuon candidates, after selection"
        histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
        
        for muTag in ["mu1","mu2"]:
            
            ARRAY=data[f'selected_dimu_{muTag}'].pt
            BINS=np.arange(0.0,40.0,0.2)
            NAME=f"dimuons_{muTag}_pt"
            DOC=f"pt {muTag} of  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
            
            ARRAY=data[f'selected_dimu_{muTag}'].eta
            BINS=np.arange(-3.0,3.0,0.1)
            NAME=f"dimuons_{muTag}_eta"
            DOC=f"eta {muTag} of  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
            
            ARRAY=data[f'selected_dimu_{muTag}'].d0
            BINS=np.arange(0.0,4.0,0.02)
            NAME=f"dimuons_{muTag}_d0"
            DOC=f"d0 / {muTag} of  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
            
            ARRAY=data[f'selected_dimu_{muTag}'].iso
            BINS=np.arange(0.0,10.0,0.2)
            NAME=f"dimuons_{muTag}_tkIsoSum"
            DOC=f"tk Iso Sum of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

            ARRAY=data[f'selected_dimu_{muTag}'].iso/data[f'selected_dimu_{muTag}'].pt
            BINS=np.arange(0.0,5.0,0.1)
            NAME=f"dimuons_{muTag}_tkRelIso"
            DOC=f"tk relative Iso  of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

            ARRAY=data[f'selected_dimu_{muTag}'].pf_iso_sum
            BINS=np.arange(0.0,10.0,0.2)
            NAME=f"dimuons_{muTag}_pfIsoSum"
            DOC=f"pf Iso Sum of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

            ARRAY=data[f'selected_dimu_{muTag}'].pf_rel_iso
            BINS=np.arange(0.0,5.0,0.1)
            NAME=f"dimuons_{muTag}_pfRelIso"
            DOC=f"pf relative Iso  of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
            
            ARRAY=data[f'selected_dimu_{muTag}'].pf_nPartInCone
            BINS=np.arange(-0.5,50.6,1.0)
            NAME=f"dimuons_{muTag}_pfNParticleInCone"
            DOC=f"number of pf candidates in cone of of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

            ARRAY=data[f'selected_dimu_{muTag}'].pfNeutral_iso_sum
            BINS=np.arange(0.0,10.0,0.2)
            NAME=f"dimuons_{muTag}_pfNeutralIsoSum"
            DOC=f"pfNeutral Iso Sum of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

            ARRAY=data[f'selected_dimu_{muTag}'].pfNeutral_rel_iso
            BINS=np.arange(0.0,5.0,0.1)
            NAME=f"dimuons_{muTag}_pfNeutralRelIso"
            DOC=f"pfNeutral relative Iso  of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
            
            ARRAY=data[f'selected_dimu_{muTag}'].pfNeutral_nPartInCone
            BINS=np.arange(-0.5,50.6,1.0)
            NAME=f"dimuons_{muTag}_pfNeutralNParticleInCone"
            DOC=f"number of pfNeutral candidates in cone of of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

            ARRAY=data[f'selected_dimu_{muTag}'].pfCharged_iso_sum
            BINS=np.arange(0.0,10.0,0.2)
            NAME=f"dimuons_{muTag}_pfChargedIsoSum"
            DOC=f"pfCharged Iso Sum of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

            ARRAY=data[f'selected_dimu_{muTag}'].pfCharged_rel_iso
            BINS=np.arange(0.0,5.0,0.1)
            NAME=f"dimuons_{muTag}_pfChargedRelIso"
            DOC=f"pfCharged relative Iso  of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
            
            ARRAY=data[f'selected_dimu_{muTag}'].pfCharged_nPartInCone
            BINS=np.arange(-0.5,50.6,1.0)
            NAME=f"dimuons_{muTag}_pfChargedNParticleInCone"
            DOC=f"number of pfCharged candidates in cone of of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
            
            ARRAY=data[f'selected_dimu_{muTag}'].puppi_iso_sum
            BINS=np.arange(0.0,10.0,0.2)
            NAME=f"dimuons_{muTag}_puppiIsoSum"
            DOC=f"puppi Iso Sum of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

            ARRAY=data[f'selected_dimu_{muTag}'].puppi_rel_iso
            BINS=np.arange(0.0,5.0,0.1)
            NAME=f"dimuons_{muTag}_puppiRelIso"
            DOC=f"puppi relative Iso  of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)
            
            ARRAY=data[f'selected_dimu_{muTag}'].puppi_nPartInCone
            BINS=np.arange(-0.5,50.6,1.0)
            NAME=f"dimuons_{muTag}_puppiNParticleInCone"
            DOC=f"number of puppi candidates in cone of of {muTag}  dimuon candidates, after selection"
            histStore[NAME]  =getHistograms(ARRAY,bins=BINS,name=NAME,doc=DOC)

        ofname=f"{prefix}/out_data_{ky}_{tag}.json"   
        os.system(f"mkdir  -p {prefix}")
        print(f"exporting to {ofname}")
    
        config_dict["HISTSTORE"]=histStore

        with open(ofname,"w") as f:
            json.dump(config_dict,f,indent=4)
        
        #if args.exportFile:
    if False:
        ofname=f"{prefix}/out_data_{ky}_{tag}.pkl"   
        dataExport={}
        for ky in categorizedData:
            data=categorizedData[ky]
            data=cleanEmptyDimuons(data,cleanEvents=True)
            dataExport[ky]=data[branchesToStore]
        print("Exporting file : ",ofname)
        with open(ofname,'wb') as f:
            pkl.dump(dataExport,f)


if __name__=='__main__':
    main()

