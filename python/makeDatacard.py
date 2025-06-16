import sys,os
import json,argparse
import  tabulate
import glob
import numpy as np
from util import getYieldInHumanReadableString

lumi=400e3
lumi=1e3

with open("data/xsdb.json") as f:
    xsData=json.load(f)

parser = argparse.ArgumentParser()
parser.add_argument("-v","--version", help="version",default='-1')
parser.add_argument("-g","--getInclusiveShape", help="use inclusive shape",default=False,action='store_true')
args = parser.parse_args()

masssPoints={}
masssPoints['DARK_PHOTON_M1']=(1.0-0.2,1.0+0.2)
masssPoints['DARK_PHOTON_M2']=(2.0-0.2,2.0+0.2)
masssPoints['DARK_PHOTON_M4']=(4.0-0.2,4.0+0.2)
masssPoints['DARK_PHOTON_M8']=(8.0-0.2,8.0+0.2)

def pad(s,length=12,sp=" ", ):
    nspace = max(1, length - len(s))
    return s + (sp * nspace)
    
def getCardHeader(input_data,imax,jmax,cats=None):
    if cats is None:
        cats=input_data['categories']

    outDCard=f"{pad('',20,'-')}\nimax {imax}\njmax {jmax}\nkmax *\n{pad('',20,'-')}\n"
    nmax=max([len(c) for c in cats])+5
    txt="".join(f"{c:>{nmax}}" for c in cats )
    outDCard+="bin         "+txt
    outDCard+="\n"
    obs=[]
    
    for c in cats:
        
        if 'observation' in input_data:
            obs.append( input_data['observation'][c]  )
        else: 
            obs.append( -1 )
            
    txt="".join([f"{ob:>{nmax}}" for ob in obs])
    outDCard+="observation "+txt
    outDCard+="\n"
    outDCard+=pad("",20,"-")+"\n"
    return outDCard

def getProcRateMap(cat_proc_rate_map,prc_idx_map,syst_car_proc_map={}):
    od=""
    l1="bin          "
    l2="process      "
    l3="process      "
    l4="rate         "
    ls="lumi    lnN  "
    lTPL="@@SYST    lnN  "
    
    nmax=0
    for cat in cat_proc_rate_map:
        nmax=max(nmax,len(cat))
        for proc in cat_proc_rate_map[cat]:
            rate=cat_proc_rate_map[cat][proc]
            nmax=max(nmax,len(proc))
            nmax=max(nmax,len(f"{rate:.1f}"))
    nmax+=5
    allSystLines={syst:f"{lTPL.replace('@@SYST',syst):>12}" for syst in syst_car_proc_map}
    for cat in cat_proc_rate_map:
        for proc in cat_proc_rate_map[cat]:
            pidx=prc_idx_map[proc]
            rate=cat_proc_rate_map[cat][proc]
            l1+= f" {cat :>{nmax}}"
            l2+= f" {proc:>{nmax}}"
            l3+= f" {pidx:>{nmax}}"
            # l4+= f" {rate:>{nmax}}"
            l4+= f" {rate:>{nmax}.1f}"
            ls+= f" {'1.1':>{nmax}}"
            for syst in  syst_car_proc_map :
                if  cat in syst_car_proc_map[syst]:
                    if proc in syst_car_proc_map[syst][cat]:
                        x=f'{syst_car_proc_map[syst][cat][proc]:.3f}'
                        allSystLines[syst]+= f" {x:>{nmax}}"
                    else:
                        allSystLines[syst]+= f" {' -':>{nmax}}"
                else:
                    allSystLines[syst]+= f" {' -':>{nmax}}"
    
    od+=l1 + "\n"
    od+=l2 + "\n"
    od+=l3 + "\n"
    od+=l4 + "\n"
    od+=pad("",20,"-")+"\n"
    od+=ls + "\n"
    for ky in allSystLines:
        od+=allSystLines[ky] + "\n"

    return od

def getDatacard(input_data,cats=None):
    if cats is None:
        cats=input_data['categories']
    imax=len(cats)
    jmax=len(input_data['backgrounds']+input_data['signal'])-1
    
    cat_proc_rate_map={}
    cat_proc_syst_map={}
    prc_idx_map={}
    for i,sig in enumerate(input_data['signal']):
        prc_idx_map[sig]=-1*i
    for i,bkg in enumerate(input_data['backgrounds']):
        prc_idx_map[bkg]=1+i
    syst_car_proc_map={}
    for prc in input_data['signal']+input_data['backgrounds']:
        for cat in cats: # input_data['rate'][prc]
            if cat not in cat_proc_rate_map:
                cat_proc_rate_map[cat]={}
            cat_proc_rate_map[cat][prc]=input_data['rate'][prc][cat]
            for syst in input_data['systematics'][prc][cat]:
                if syst not in syst_car_proc_map:
                    syst_car_proc_map[syst]={}
                if cat not in syst_car_proc_map[syst]:
                    syst_car_proc_map[syst][cat]={}
                syst_car_proc_map[syst][cat][prc]=input_data['systematics'][prc][cat][syst]
            
                
    
    outDCard=getCardHeader(input_data,imax,jmax,cats)
    od=getProcRateMap(cat_proc_rate_map,prc_idx_map,syst_car_proc_map)

    dcard_text=outDCard+od
    
    return dcard_text

ver=args.version
base  =f'results/analysis/{ver}/'
prefix=f'results/datacards/{ver}/'
if args.getInclusiveShape:
    prefix=f'results/datacards/{ver}/WithInclusiveShape/'

print("Looking at files in ",base)
print("Exporting to ",prefix)

procFolders=glob.glob(base+'/*')
procs=[fd.split("/")[-1] for fd in procFolders]
prc=procs[0]
jFiles=glob.glob(base+'/'+procs[0]+"/*.json")
cats=[jl.split('/')[-1].replace('out_data_','').replace(f'_{prc}','').replace('.json','') for jl in jFiles ]
allCats=list(cats)
cats=[ct for ct in cats if ct not in ['inclusive'] ]
print(cats)

data_dict_catWise={}
data_dict_procWise={}
for cat in allCats:
    print("= "*10," Processing category ",cat)
    for proc in procs:
        srchstr=f"{base}/{proc}/*_{cat}_{proc}*.json"
        flist=glob.glob(srchstr)
        if len(flist)!=1:
            print(f"Ambiguty in the file search ! search string : {srchstr}")
            exit(1)
        fl=flist[0]
        with open(fl) as f:
            summary=json.load(f)
        if cat not in data_dict_catWise:
            data_dict_catWise[cat]={}
        if proc not in data_dict_procWise:
            data_dict_procWise[proc]={}
        data_dict_catWise[cat][proc] =summary
        data_dict_procWise[proc][cat]=summary
    print()

KHY='dimuonMassFineBin50Mev'
signalProcs=['DARK_PHOTON_M1','DARK_PHOTON_M2','DARK_PHOTON_M4','DARK_PHOTON_M8']
bkgProcs=['DY1Jet_Mll_1To10','MinBias']

allDaracardData={}
failedModelData=[]
for sproc in signalProcs:
    input_data={}
    input_data['categories']=allCats
    input_data['backgrounds']=bkgProcs
    input_data['signal']=[sproc]
    xmin,xmax=masssPoints[sproc]
    print(f"For {sproc} : [{xmin},{xmax}] GeV window")
    rateMap={}
    systMap={}
    nBkgExpected=0.0
    for proc in bkgProcs+[sproc]:
        xs=xsData[proc]
        rateMap[proc]={}
        systMap[proc]={}
        for cat in data_dict_procWise[proc]:
            histStore=data_dict_procWise[proc][cat]
            ntot=histStore['NEVENTS']
            binx=np.array( histStore["HISTSTORE"][KHY]['bins'] )
            mask=np.logical_and( binx < xmax , binx >= xmin)
            counts=np.array(histStore["HISTSTORE"][KHY]['counts'])
            npass=np.sum(counts[mask[:-1]])
            nExpected     = npass/ntot * lumi * xs
            nExpected_err = nExpected/np.sqrt(npass+1e-15)

            if args.getInclusiveShape:
                histStore_=histStore
                histStore=data_dict_procWise[proc]['inclusive']
                binx=np.array( histStore["HISTSTORE"][KHY]['bins'] )
                mask=np.logical_and( binx < xmax , binx >= xmin)
                counts=np.array(histStore["HISTSTORE"][KHY]['counts'])
                npass=np.sum(counts[mask[:-1]])
                eff = npass/ntot
                
                histStore=histStore_
                counts     = np.array(histStore["HISTSTORE"][KHY]['counts'])
                mask       = np.logical_and( binx < 1e6 , binx >0)
                npass      = np.sum(counts[mask[:-1]])
                nExpected  = npass/ntot * lumi * xs

                nExpected_err = nExpected/np.sqrt(npass+1e-15)


            rateMap[proc][cat]=nExpected
            systMap[proc][cat]={}
            systMap[proc][cat]['norm']=1.0+nExpected_err/(nExpected+1e-16)

            #print(cat,proc,npass,nExpected)
    catsToTake=[]
    for cat in input_data['categories']:
        nBkgExpected=0.0
        for proc in bkgProcs:
            nBkgExpected+=rateMap[proc][cat]
        sig=rateMap[sproc][cat]/np.sqrt(nBkgExpected)
        s=getYieldInHumanReadableString(rateMap[sproc][cat])
        b=getYieldInHumanReadableString(nBkgExpected)
        print("  -> ",sproc,",",cat,f"  :  expected significnce :  {sig:.4f}  |   S ={s} / B ={b}" )
        if nBkgExpected==0.0 :
            print(f"For proc {sproc} in {cat} , background modeling failed ! ")
            failedModelData.append( [ cat, sproc, 'bkg modeling failed' ] )
            continue
        if rateMap[sproc][cat]==0.0 :
            print(f"For proc {sproc} in {cat} , signal modeling failed ! ")
            failedModelData.append( [ cat, sproc, 'sig modeling failed' ] )
            continue
        catsToTake.append(cat)
        

    input_data['categories']=catsToTake
    input_data['rate']=rateMap
    input_data['systematics']=systMap
    allDaracardData[sproc]=input_data 

cmd=f'mkdir -p {prefix}/perCat/'
os.system(cmd)
for ky in allDaracardData:
    print("\n"+ky)
    dcTxt=getDatacard(allDaracardData[ky])
    ofname=f"{prefix}/{ky}.txt" 
    print("Exporting : ",ofname)
    with open(ofname,'w') as f:
        f.write(dcTxt)
    k=0
    for i,cat in  enumerate(allDaracardData[ky]['categories']):
        dcTxt=getDatacard(allDaracardData[ky],[cat])
        ofname=f"{prefix}/perCat/{ky}_{cat}.txt" 
        print(f"\r[{i+1}] Exporting : ",ofname,allDaracardData[ky]['rate'][ky][cat], "    ",end="")
        with open(ofname,'w') as f:
            f.write(dcTxt)
    print()
for itm in failedModelData:
    print(itm)
