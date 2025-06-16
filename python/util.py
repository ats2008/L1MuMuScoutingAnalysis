import uproot as urt
import json
import numpy as np
import awkward as ak
import vector
vector.register_awkward()

MASS_MU=0.1056583755

signalSamples=['DARK_PHOTON_M1', 'DARK_PHOTON_M2', 'DARK_PHOTON_M4', 'DARK_PHOTON_M8']
bkgSamples=['DY_Mll_10To50','DY1Jet_Mll_1To10',  'MinBias'] # 'QCD15To20_Mu5Enriched']

labelMap={}
labelMap['MinBias']="Min-Bias"
labelMap['DARK_PHOTON_M1']="DPh(M=1)"
labelMap['DARK_PHOTON_M2']="DPh(M=2)"
labelMap['DARK_PHOTON_M4']="DPh(M=4)"
labelMap['DARK_PHOTON_M8']="DPh(M=8)"
labelMap['DY1Jet_Mll_1To10']="DY( m$_{ll}\in$1-10 )"
labelMap['DY_Mll_10To50']="DY( m$_{ll}\in$10-50 )"

def getLabelFromTag(tag):
    if tag in labelMap:
        return labelMap[tag]
    return tag
def getYieldInHumanReadableString(nExpected,nExpected_err=None):
    if nExpected > 1e8:
        yld=f"{nExpected/1e9:.2f} B"
        if nExpected_err is not None:
            yld=yld[:-1]+f"+/- {nExpected_err/1e9:.2f} B"
    elif nExpected > 1e5:
        yld=f"{nExpected/1e6:.2f} M"
        if nExpected_err is not None:
            yld=yld[:-1]+f"+/- {nExpected_err/1e6:.2f} M"
    elif nExpected > 1e3:
        yld=f"{nExpected/1e3:.2f} K"
        if nExpected_err is not None:
            yld=yld[:-1]+f"+/- {nExpected_err/1e9:.2f} K"
    else :
        yld=f"{nExpected:.2f} "
        if nExpected_err is not None:
            yld=yld[:-1]+f"+/- {nExpected_err:.2f} "
    return yld



def getHistograms(array,bins,name="hist",doc="",isFaltArray=False):
    if not isFaltArray:
        if array.ndim > 1:
            array=ak.ravel(array)
    hist={}
    c,b=np.histogram(array,bins=bins)
    hist["NAME"]=name
    hist["DOC"]=doc
    hist["counts"]=[ float(i) for i in c]
    hist["bins"]  =[ float(i) for i in b]
    return hist

def addGenMatchedMuonInformation(dataStore,pTMinMuon=2.0,etaMax=2.4):
    print("Adding GEN-MATCH infomation ( GEN vs TkMuon ) ")
    for tag in dataStore:
        print("procesing ",tag)
        data=dataStore[tag]
        mask= (np.abs(data.GenCands_pdgId) == 13) & (np.abs(data.GenCands_eta) < 2.4) & (data.GenCands_pt > pTMinMuon)
        gen_muons= ak.zip(
            {
                'pt' : data.GenCands_pt[mask],
                'eta' : data.GenCands_eta[mask],
                'phi' : data.GenCands_phi[mask],
                'mass' : data.GenCands_mass[mask]
            },
             with_name="Momentum4D"
        )

        # adding mass var for tkMuon 4 Momentum
        data['L1gmtTkMuon_mass']=data['L1gmtTkMuon_pt']*0.0 +MASS_MU
        tkMuons=ak.zip( {ky.split("_")[1] : data[ky] for ky in data.fields  if 'L1gmtTkMuon_' in ky} ,  with_name="Momentum4D")
        mask=ak.num(tkMuons) <1
        
        # making proxy 0 momentum muons for drMatch  //  NOT NEEDED
        #zero_fourvec=ak.zip({'pt':0,'eta':0.0,'phi':0.0,'mass':0.0},with_name='Momentum4D')
        #for ky in tkMuons.fields:
        #    zero_fourvec[ky]=0.0
        #tkMuons_=ak.fill_none(ak.pad_none(tkMuons,1,axis=-1),zero_fourvec)
        tkMuons_=tkMuons
        # Pairing
        gm,tm=ak.unzip(ak.cartesian((gen_muons, tkMuons_),nested=True))
        
        # matching
        dR=gm.deltaR(tm)
        midx=ak.argmin(dR,axis=-1,keepdims=True)
        drMins=dR[midx]
        matchedMask=drMins<0.1
        
        l1SignalMuons=tm[midx][matchedMask]
        l1SignalMuons=ak.flatten(l1SignalMuons,axis=-1)
    
        cmb=ak.combinations( l1SignalMuons, 2, axis=-1, replacement=False)
        mu1_dimu_,mu2_dimu_=cmb['0'],cmb['1']
        
        # Assignment
        dimu=mu1_dimu_+mu2_dimu_
        data['matchedDimuons']=dimu
    
        data['matchedDimuons_mu1']=mu1_dimu_
        data['matchedDimuons_mu2']=mu2_dimu_
    return dataStore   
 
def addAllDimuonCombinations(dataStore):
    print("Making all dimuon candidates !")
    for tag in dataStore: 
        print("   > Processing ",tag)
        data=dataStore[tag]
    
        data['L1gmtTkMuon_mass']=data['L1gmtTkMuon_pt']*0.0 +MASS_MU
        tkMuons=ak.zip( {ky.split("_")[1] : data[ky] for ky in data.fields  if 'L1gmtTkMuon_' in ky} ,  with_name="Momentum4D")
    
        #zero_fourvec=ak.zip({'pt':0,'eta':0.0,'phi':0.0,'mass':0.0},with_name='Momentum4D')
        #for ky in tkMuons.fields:
        #    zero_fourvec[ky]=0.0
        #tkMuons_=ak.fill_none(ak.pad_none(tkMuons,1,axis=-1),zero_fourvec)
        tkMuons_=tkMuons
    
        dimu_=ak.combinations(tkMuons_, 2)
        ch_prod=dimu_['0'].charge*dimu_['1'].charge
        deltaZMask = np.abs(dimu_['0'].z0 - dimu_['1'].z0) < 0.5
        msk= (ch_prod < 0) & deltaZMask
        dimu_=dimu_[msk]
        dimu_m1=dimu_['0']
        dimu_m2=dimu_['1']
        dimu=dimu_m1 + dimu_m2
        ptd=dimu.pt
        srtIdx=ak.argsort(ptd*-1)
        data["all_dimuons"]=dimu[srtIdx]
        data["all_dimuons_m1"]=dimu_m1[srtIdx]
        data["all_dimuons_m2"]=dimu_m2[srtIdx]

    return dataStore
    
def addSelectedDimuonsWithIsoCut(dataStore,pTMinMuon=0.0,etaMax=3.0,addTkIsoMask=True,puppiIsoThr=None):
    print("Making the selected dimuon pairs ")
    for tag in dataStore:
        print(" Processing  ",tag)
        data = dataStore[tag]
        if puppiIsoThr is not None:
            mu1_iso_mask = data.all_dimuons_m1.puppi_rel_iso < puppiIsoThr
            mu2_iso_mask = data.all_dimuons_m2.puppi_rel_iso < puppiIsoThr
            mask = mask & mu1_iso_mask
            mask = mask & mu2_iso_mask


def addSelectedDimuons(dataStore,pTMinMuon=0.0,etaMax=3.0,addTkIsoMask=True,doDRSelection=False,d0Max=-1e4,puppiIsoThr=None,pfIsoThr=None):
    dataStore=addAllDimuonCombinations(dataStore)

    print("Making the selected dimuon pairs ")
    for tag in dataStore:
        print(" Processing  ",tag)
        data = dataStore[tag]
    
        pt1_overM_mask = data.all_dimuons_m1.pt / data.all_dimuons.mass > 0.25
        pt2_overM_mask = data.all_dimuons_m2.pt / data.all_dimuons.mass > 0.25
    
        mu1_id_mask = data.all_dimuons_m1.tightId > 0.5
        mu2_id_mask = data.all_dimuons_m2.tightId > 0.5

        mu1_pt_mask = data.all_dimuons_m1.pt > pTMinMuon
        mu2_pt_mask = data.all_dimuons_m2.pt > pTMinMuon

        mu1_eta_mask = np.abs(data.all_dimuons_m1.eta) < etaMax
        mu2_eta_mask = np.abs(data.all_dimuons_m2.eta) < etaMax

        deltaZMask = np.abs(data.all_dimuons_m1.z0 - data.all_dimuons_m2.z0) < 0.5
    
        mask = pt1_overM_mask & pt2_overM_mask
        mask = mask & mu1_id_mask
        mask = mask & mu2_id_mask
        
        mask = mask & mu1_pt_mask
        mask = mask & mu2_pt_mask
        
        mask = mask & mu1_eta_mask
        mask = mask & mu2_eta_mask

        mask = mask & deltaZMask
        
        if addTkIsoMask:
            mu1_iso_mask = data.all_dimuons_m1.iso < 0.5
            mu2_iso_mask = data.all_dimuons_m2.iso < 0.5
            mask = mask & mu1_iso_mask
            mask = mask & mu2_iso_mask
 
        if d0Max > 0:
            mu1_d0_mask = data.all_dimuons_m1.d0 < d0Max
            mu2_d0_mask = data.all_dimuons_m2.d0 < d0Max
            mask = mask & mu1_d0_mask
            mask = mask & mu2_d0_mask
            
        if doDRSelection:
            mumu_deltaR_mask = data.all_dimuons_m1.deltaR( data.all_dimuons_m2  ) <  1.3
            mask = mask & mumu_deltaR_mask
    
        data["selected_dimu"]    = data["all_dimuons"][mask]
        data["selected_dimu_mu1"]= data["all_dimuons_m1"][mask]
        data["selected_dimu_mu2"]= data["all_dimuons_m2"][mask]

        if puppiIsoThr is not None:
            if 'puppi_rel_iso' not in data.selected_dimu_mu1.fields:
                data=addPuppiIsolation({'sel':data})['sel']
            mu1_iso_mask = data.selected_dimu_mu1.puppi_rel_iso < puppiIsoThr
            mu2_iso_mask = data.selected_dimu_mu2.puppi_rel_iso < puppiIsoThr
            mask =  mu1_iso_mask & mu2_iso_mask
            data["selected_dimu"]    =     data["selected_dimu"][mask]     
            data["selected_dimu_mu1"]=     data["selected_dimu_mu1"][mask] 
            data["selected_dimu_mu2"]=     data["selected_dimu_mu2"][mask] 

        if pfIsoThr is not None:
            if 'pf_rel_iso' not in data.selected_dimu_mu1.fields:
                data=addPFIsolations({'sel':data})['sel']
            mu1_iso_mask = data.selected_dimu_mu1.pf_rel_iso < pfIsoThr
            mu2_iso_mask = data.selected_dimu_mu2.pf_rel_iso < pfIsoThr
            mask =  mu1_iso_mask & mu2_iso_mask
            data["selected_dimu"]    =     data["selected_dimu"][mask]     
            data["selected_dimu_mu1"]=     data["selected_dimu_mu1"][mask] 
            data["selected_dimu_mu2"]=     data["selected_dimu_mu2"][mask] 



        dimu=data["selected_dimu"]
        ptd=dimu.pt
        srtIdx=ak.argsort(ptd*-1)
        data["selected_dimu"]   =data["selected_dimu"][srtIdx]
        data["selected_dimu_mu1"]=data["selected_dimu_mu1"][srtIdx]
        data["selected_dimu_mu2"]=data["selected_dimu_mu2"][srtIdx]

        for col in ["selected_dimu","selected_dimu_mu1","selected_dimu_mu2"]:
            zero_fourvec=ak.zip({'pt':0,'eta':1e3,'phi':1e3,'mass':0.0},with_name='Momentum4D')
            for ky in data[col].fields:
                if ky not in ['eta','phi']:
                    zero_fourvec[ky]=0.0
            data[col]=ak.fill_none(ak.pad_none(data[col],1,axis=-1),zero_fourvec)

        data["selected_dimu"]    =ak.unflatten(data["selected_dimu"][:,0],counts=1)
        data["selected_dimu_mu1"]=ak.unflatten(data["selected_dimu_mu1"][:,0],counts=1)
        data["selected_dimu_mu2"]=ak.unflatten(data["selected_dimu_mu2"][:,0],counts=1)

    return dataStore

    
def addSelectedDimuonsForJpsi(dataStore):
    dataStore=addAllDimuonCombinations(dataStore)

    print("Making the selected dimuon pairs for Jpsi ")
    for tag in dataStore:
        print(" Processing  ",tag)
        data = dataStore[tag]
    
        mu1_id_mask = data.all_dimuons_m1.tightId > 0.5
        mu2_id_mask = data.all_dimuons_m2.tightId > 0.5
        
        mumu_deltaR_mask = data.all_dimuons_m1.deltaR( data.all_dimuons_m2  ) <  1.5

        mu1_iso_mask = data.all_dimuons_m1.iso < 0.5
        mu2_iso_mask = data.all_dimuons_m2.iso < 0.5
        
        deltaZMask = np.abs(data.all_dimuons_m1.z0 - data.all_dimuons_m2.z0) < 0.5
        
        mask =  data.all_dimuons_m1.tightId > -1e4
        mask = mask & mu1_id_mask
        mask = mask & mu2_id_mask
        #mask = mask & mu1_iso_mask
        #mask = mask & mu2_iso_mask
        mask = mask & deltaZMask
        mask = mask & mumu_deltaR_mask
    
        data["selected_dimu"]= data["all_dimuons"][mask]
        data["selected_dimu_mu1"]= data["all_dimuons_m1"][mask]
        data["selected_dimu_mu2"]= data["all_dimuons_m2"][mask]
        for col in ["selected_dimu","selected_dimu_mu1","selected_dimu_mu2"]:
            zero_fourvec=ak.zip({'pt':0,'eta':1e3,'phi':1e3,'mass':0.0},with_name='Momentum4D')
            for ky in data[col].fields:
                if ky not in ['eta','phi']:
                    zero_fourvec[ky]=0.0
            data[col]=ak.fill_none(ak.pad_none(data[col],1,axis=-1),zero_fourvec)
    
    return dataStore

#def addPFIsolationVariables(dataStore,muBaseString="matchedDimuons"):
#    for tag in dataStore:
#        data=dataStore[tag]
#        for muBase in [muBaseString]:
#            print("Processing ",muBase)
#            if muBase not in data.fields:
#                print("      Skipping ",muBase," for ",tag)
#                continue
#            
#            pCands=ak.zip( {ky.split("_")[1] : data[ky] for ky in data.fields  if 'L1PFCands_' in ky} ,
#                          with_name="Momentum4D")
#            
#            for muTag,muOtag in zip(["mu1","mu2"],["mu2","mu1"]):
#                print(f"        > Doing {muBase}_{muTag}",)
#                muX=data[f'{muBase}_{muTag}']
#                muO=data[f'{muBase}_{muOTag}']
#        #         zero_fourvec=ak.zip({'pt':0,'eta':0.0,'phi':0.0,'mass':0.0},with_name='Momentum4D')
#        #         for ky in pfCands.fields:
#        #             zero_fourvec[ky]=0.0
#        #         pCands=ak.fill_none(ak.pad_none(pCands,1,axis=-1),zero_fourvec)
#                mu ,pc =ak.unzip(ak.cartesian((muX, pCands),nested=True))
#                muo,pco=ak.unzip(ak.cartesian((muO, pCands),nested=True))
#    
#                dR=mu.deltaR(pc)
#                dRO=muo.deltaR(pco)
#                dr_mask = (dR > 0.05) & ( dR < 0.3 ) & ( dRO > 0.02 )
#
#                delataZ = np.abs(pc.z0 - mu.z0)
#                ch_selection_mask = (np.abs(pc.charge) > 0.5 ) & (delataZ < 0.5)
#                nu_selection_mask = (np.abs(pc.charge) < 0.5 )
#                selection_mask_in_cone = ( nu_selection_mask | ch_selection_mask ) & dr_mask
#    
#                pc_inCone = pc[selection_mask_in_cone] 
#                nPInCone=ak.num(pc_inCone,axis=-1)
#    
#                iso_contibution = pc_inCone.pt #*pc_inCone.puppiWeight
#                iso_sum = ak.sum(iso_contibution,axis=-1)
#                #mu_poi= mu[...,0]
#                mu_poi= muX
#                rel_iso=iso_sum/(mu_poi.pt +1e-10)
#    #            print(rel_iso)
#    
#                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],nPInCone,"pf_nPartInCone")
#                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],iso_sum ,"pf_iso_sum")
#                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],rel_iso ,"pf_rel_iso")
#
#    return dataStore



#def addPUPPIIsolationVariables(dataStore,muBaseString="matchedDimuons"):
#    for tag in dataStore:
#        data=dataStore[tag]
#        
#        for muBase in ['matchedDimuons','selected_dimu']:
#            print("Processing ",muBase)
#            if muBase not in data.fields:
#                print("      Skipping ",muBase," for ",tag)
#                continue
#            
#            pCands=ak.zip( {ky.split("_")[1] : data[ky] for ky in data.fields  if 'L1PuppiCands_' in ky} ,
#                          with_name="Momentum4D")
#            for muTag,muOtag in zip(["mu1","mu2"],["mu2","mu1"]):
#                print(f"        > Doing {muBase}_{muTag}",)
#                muX=data[f'{muBase}_{muTag}']
#                muO=data[f'{muBase}_{muOTag}']
#                if 'puppi_rel_iso' in muX.fields:
#                    continue
#        #         zero_fourvec=ak.zip({'pt':0,'eta':0.0,'phi':0.0,'mass':0.0},with_name='Momentum4D')
#        #         for ky in pfCands.fields:
#        #             zero_fourvec[ky]=0.0
#        #         pCands=ak.fill_none(ak.pad_none(pCands,1,axis=-1),zero_fourvec)
#                mu ,pc =ak.unzip(ak.cartesian((muX, pCands),nested=True))
#                muo,pco=ak.unzip(ak.cartesian((muO, pCands),nested=True))
#    
#                dR=mu.deltaR(pc)
#                dRO=muo.deltaR(pco)
#                dr_mask = (dR > 0.05) & ( dR < 0.3 ) & ( dRO > 0.02 )
#
#                pc_inCone = pc[dr_mask] 
#                nPInCone=ak.num(pc_inCone,axis=-1)
#    
#                iso_contibution = pc_inCone.pt *pc_inCone.puppiWeight
#                iso_sum = ak.sum(iso_contibution,axis=-1)
#                #mu_poi= mu[...,0]
#                mu_poi= muX
#                rel_iso=iso_sum/(mu_poi.pt+1e-9)
#    
#                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],nPInCone,"puppi_nPartInCone")
#                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],iso_sum,"puppi_iso_sum")
#                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],rel_iso,"puppi_rel_iso")
#    
#    return dataStore

def getCategories(data,cfg):
#        data["selected_dimu"]= data["all_dimuons"][mask]
#        data["selected_dimu_mu1"]= data["all_dimuons_m1"][mask]
#        data["selected_dimu_mu2"]= data["all_dimuons_m2"][mask]
    
    catMask=data["selected_dimu_mu1"].pt > -1
    
    print(f"All {ak.sum(catMask)}")
    if 'min_pt_mu1' in cfg:
        ptS=cfg['min_pt_mu1']
        catMask = catMask & ( data["selected_dimu_mu1"].pt >= ptS )
        print(f"min_pt_mu1 {ak.sum(catMask)}  {ptS}")

    if 'min_pt_mu2' in cfg:
        ptS=cfg['min_pt_mu2']
        catMask = catMask & ( data["selected_dimu_mu2"].pt >= ptS )
        print(f"min_pt_mu2 {ak.sum(catMask)} {ptS}")

    if 'min_eta_mu1' in cfg:
        etaS=cfg['min_eta_mu1']
        catMask = catMask & ( np.abs(data["selected_dimu_mu1"].eta) >= etaS )
        print(f"min_eta_mu1 {ak.sum(catMask)} {etaS}")

    if 'min_eta_mu2' in cfg:
        etaS=cfg['min_eta_mu2']
        catMask = catMask & ( np.abs(data["selected_dimu_mu2"].eta) >= etaS )
        print(f"min_eta_mu2 {ak.sum(catMask)}  {etaS}")

    if 'max_eta_mu1' in cfg:
        etaS=cfg['max_eta_mu1']
        catMask = catMask & ( np.abs(data["selected_dimu_mu1"].eta) < etaS )
        print(f"max_eta_mu1 {ak.sum(catMask)} {etaS}")

    if 'max_eta_mu2' in cfg:
        etaS=cfg['max_eta_mu2']
        catMask = catMask & ( np.abs(data["selected_dimu_mu2"].eta) < etaS )
        print(f"max_eta_mu2 {ak.sum(catMask)}  {etaS}")

    if 'min_dimu_pt' in cfg:
        ptS=cfg['min_dimu_pt']
        catMask = catMask & ( data["selected_dimu"].pt >ptS )
        print(f"min_dimu_pt {ak.sum(catMask)}  {ptS}")

    if 'min_dimu_dr' in cfg:
        drMin=cfg['min_dimu_dr']
        dr=data["selected_dimu_mu1"].deltaR( data["selected_dimu_mu2"]  )
        catMask = catMask & ( dr >= drMin )
    
    if 'max_dimu_dr' in cfg:
        drMax=cfg['max_dimu_dr']
        dr=data["selected_dimu_mu1"].deltaR( data["selected_dimu_mu2"]  )
        catMask = catMask & ( dr < drMax )
    
    dimu_cat_event_mask =  ak.sum(catMask , axis = -1) > 0
    dimu_cat_event_inveted_mask = np.logical_not( dimu_cat_event_mask  )
    print(ak.sum(dimu_cat_event_mask),ak.sum(dimu_cat_event_inveted_mask))
    return data[dimu_cat_event_mask],data[dimu_cat_event_inveted_mask]

def cleanEmptyDimuons(data,cleanEvents=False):    
    if 'all_dimuons' in data.fields:
        mask = data['all_dimuons'].pt > 0.0
        data['all_dimuons']=data['all_dimuons'][mask]
        data['all_dimuons_m1']=data['all_dimuons_m1'][mask]
        data['all_dimuons_m2']=data['all_dimuons_m2'][mask]
    
    if 'matched_dimu' in data.fields:
        mask = ( data['matched_dimu_mu1'].pt > 0.0 ) & (  data['matched_dimu_mu2'].pt > 0.0 )
        data['matched_dimu']=data['matched_dimu'][mask]
        data['matched_dimu_mu1']=data['matched_dimu_mu1'][mask]
        data['matched_dimu_mu2']=data['matched_dimu_mu2'][mask]

    if 'selected_dimu' in data.fields:
        mask = (data['selected_dimu_mu1'].pt > 0.0) & (data['selected_dimu_mu2'].pt > 0.0)
        #print(len(mask))
        #print(np.sum(mask))
        #data=data[mask]
        data['selected_dimu']    =data['selected_dimu'][mask]
        data['selected_dimu_mu1']=data['selected_dimu_mu1'][mask]
        data['selected_dimu_mu2']=data['selected_dimu_mu2'][mask]
    
    if cleanEvents:
        emask=ak.num(data['selected_dimu']) > 0
        data=data[emask]

    return data

def addPFIsolations(dataStore):
    for tag in dataStore:
        data=dataStore[tag]
        for muBase in ['matchedDimuons','selected_dimu']:
            if muBase not in data.fields:
                #print("      Skipping ",muBase," for ",tag)
                continue
            print("  Adding PF Isolation to ",muBase)
            pCands=ak.zip( {ky.split("_")[1] : data[ky] for ky in data.fields  if 'L1PFCands_' in ky} ,with_name="Momentum4D")
            for muTag,muOtag in zip(["mu1","mu2"],["mu2","mu1"]):
                print(f"        > Doing {muBase}_{muTag}",)
                muX=data[f'{muBase}_{muTag}']
                muO=data[f'{muBase}_{muOtag}']
                if 'pf_rel_iso' in muX.fields:
                    continue
                #zero_fourvec=ak.zip({'pt':0,'eta':0.0,'phi':0.0,'mass':0.0},with_name='Momentum4D')
                #for ky in pCands.fields:
                #    zero_fourvec[ky]=0.0
                #pCands=ak.fill_none(ak.pad_none(pCands,1,axis=-1),zero_fourvec)
                mu ,pc =ak.unzip(ak.cartesian((muX, pCands),nested=True))
                muo,pco=ak.unzip(ak.cartesian((muO, pCands),nested=True))
    
                dR=mu.deltaR(pc)
                dRO=muo.deltaR(pco)
                dr_mask = (dR > 0.05) & ( dR < 0.3 ) & ( dRO > 0.02 )
                
                delataZ = np.abs(pc.z0 - mu.z0)
                ch_selection_mask = (np.abs(pc.charge) > 0.5 ) & (delataZ < 0.5) & dr_mask
                nu_selection_mask = (np.abs(pc.charge) < 0.5 ) & dr_mask
                selection_mask_in_cone =  nu_selection_mask | ch_selection_mask 
                
                pc_inCone = pc[selection_mask_in_cone] 
                nPInCone=ak.num(pc_inCone,axis=-1)
                iso_contibution = pc_inCone.pt #*pc_inCone.puppiWeight
                iso_sum = ak.sum(iso_contibution,axis=-1)
                #mu_poi= mu[...,0]
                mu_poi= muX
                rel_iso=iso_sum/mu_poi.pt
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],nPInCone,"pf_nPartInCone")
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],iso_sum,"pf_iso_sum")
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],rel_iso,"pf_rel_iso")
                
                pc_inCone = pc[nu_selection_mask] 
                nPInCone=ak.num(pc_inCone,axis=-1)
                iso_contibution = pc_inCone.pt #*pc_inCone.puppiWeight
                iso_sum = ak.sum(iso_contibution,axis=-1)
                #mu_poi= mu[...,0]
                rel_iso=iso_sum/mu_poi.pt
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],nPInCone,"pfNeutral_nPartInCone")
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],iso_sum,"pfNeutral_iso_sum")
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],rel_iso,"pfNeutral_rel_iso")
 
                pc_inCone = pc[ch_selection_mask] 
                nPInCone=ak.num(pc_inCone,axis=-1)
                iso_contibution = pc_inCone.pt #*pc_inCone.puppiWeight
                iso_sum = ak.sum(iso_contibution,axis=-1)
                #mu_poi= mu[...,0]
                mu_poi= muX
                rel_iso=iso_sum/mu_poi.pt
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],nPInCone,"pfCharged_nPartInCone")
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],iso_sum,"pfCharged_iso_sum")
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],rel_iso,"pfCharged_rel_iso")
    
    return dataStore


    
def addPuppiIsolation(dataStore):
    for tag in dataStore:
        data=dataStore[tag]
        for muBase in ['matchedDimuons','selected_dimu']:
            print("Adding puppi iso to ",muBase)
            if muBase not in data.fields:
                print("      Skipping ",muBase," for ",tag)
                continue
            pCands=ak.zip( {ky.split("_")[1] : data[ky] for ky in data.fields  if 'L1PuppiCands_' in ky} ,
                          with_name="Momentum4D")
            
            for muTag,muOtag in zip(["mu1","mu2"],["mu2","mu1"]):
                print(f"        > Doing {muBase}_{muTag}",)
                muX=data[f'{muBase}_{muTag}']
                muO=data[f'{muBase}_{muOtag}']
                if 'puppi_rel_iso' in muX.fields:
                    continue
                mu ,pc =ak.unzip(ak.cartesian((muX, pCands),nested=True))
                muo,pco=ak.unzip(ak.cartesian((muO, pCands),nested=True))
    
                dR=mu.deltaR(pc)
                dRO=muo.deltaR(pco)
                dr_mask = (dR > 0.05) & ( dR < 0.3 ) & ( dRO > 0.02 )

                pc_inCone = pc[dr_mask] 
                nPInCone=ak.num(pc_inCone,axis=-1)
                iso_contibution = pc_inCone.pt *pc_inCone.puppiWeight
                iso_sum = ak.sum(iso_contibution,axis=-1)
                #mu_poi= mu[...,0]
                mu_poi= muX
                rel_iso=iso_sum/mu_poi.pt
    
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],nPInCone,"puppi_nPartInCone")
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],iso_sum,"puppi_iso_sum")
                data[f'{muBase}_{muTag}']=ak.with_field(data[f'{muBase}_{muTag}'],rel_iso,"puppi_rel_iso")
    
    return dataStore
        

