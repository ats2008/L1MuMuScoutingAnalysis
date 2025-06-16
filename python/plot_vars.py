import sys,os
import json,argparse
import  tabulate
import glob
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import util as utl

from cycler import cycler
hep.set_style('CMS')
cols10=['#3f90da','#ffa90e','#bd1f01','#94a4a2','#832db6','#a96b59','#e76300','#b9ac70','#717581','#92dadd']
default_cycler = (cycler(color=cols10))
plt.rc('axes', prop_cycle=default_cycler)

def getTH1FromNumpHist(nph,nph_w2=None):
    hist = ROOT.TH1D("","",len(nph[1])-1,np.array(nph[1],dtype=float))
    for i,freq in enumerate(nph[0]):
        hist.SetBinContent(i+1,freq)
    if nph_w2!=None:
        for i,freq_w2 in enumerate(nph_w2[0]):
            hist.SetBinError(i+1,np.sqrt(freq_w2))
    return hist

with open("data/xsdb.json") as f:
    xsData=json.load(f)


parser = argparse.ArgumentParser()
parser.add_argument("-v","--version", help="version",default='-1')
parser.add_argument("-p","--plotIdx", help="Plot index",type=int,default=-1)
args = parser.parse_args()

plotIdx=args.plotIdx

ALL_CATS=['highptMuMuCentral', 'highptMuMuForward', 'MuMuCentral', 'MuMuForeward','rest','inclusive']
CATLIST=['highptMuMuCentral','highptMuMuForward','MuMuCentral','MuMuForeward','inclusive','rest']
signalSamples=['DARK_PHOTON_M1', 'DARK_PHOTON_M2', 'DARK_PHOTON_M4', 'DARK_PHOTON_M8']
bkgSamples=['DY_Mll_10To50','DY1Jet_Mll_1To10',  'MinBias'] # 'QCD15To20_Mu5Enriched']

base='results/analysis/v2p1/'
prefix='results/plots/v2p1/'
base='results/analysis/v4p1/'
prefix='results/plots/v4p1/'
base='results/analysis/v3/'
prefix='results/plots/v3/'
base=f'results/analysis/{args.version}'
prefix=f'results/plots/{args.version}'
cats=CATLIST
os.system(f'mkdir -p {prefix}')

cmdExtras=" "

procFolders=glob.glob(base+'/*')
procs=[fd.split("/")[-1] for fd in procFolders]
data_dict_catWise={}
data_dict_procWise={}
print("Loading from ",base)
print("Procs got : ",procs)
for cat in cats:
    print("= "*10," Processing category ",cat)
    for proc in procs:
        srchstr=f"{base}/{proc}/*_{cat}_*"
        flist=glob.glob(srchstr)
        if len(flist)!=1:
            print(f"Ambiguty in the file search ! search string : {srchstr}")
            exit(1)
        fl=flist[0]
#         print(fl)
        with open(fl) as f:
            summary=json.load(f)
        if cat not in data_dict_catWise:
            data_dict_catWise[cat]={}
        if proc not in data_dict_procWise:
            data_dict_procWise[proc]={}
        data_dict_catWise[cat][proc] =summary
        data_dict_procWise[proc][cat]=summary
    print()


def getHistograms(histStore,khy,bkgSamples,signalSamples):
    histsBkg={'hist':[],'w2':[],'label':[],'bins':None}
    histsSig={'hist':[],'w2':[],'label':[],'bins':None}

    for ky in bkgSamples:
        bins=np.array(histStore[ky]['HISTSTORE'][khy]['bins'])
        counts=np.array(histStore[ky]['HISTSTORE'][khy]['counts'])
        xs=xsData[ky]
        norm=1.0*xs*400*1e3/histStore[ky]["NEVENTS"]
        w=counts*norm
        w2=w*w
#         print(histStore[ky]["NEVENTS"],ky,np.sum(w))
        histsBkg['hist'].append(w)
        histsBkg['w2'].append(w2)
        histsBkg['label'].append(utl.getLabelFromTag(ky))

    for ky in signalSamples:
        bins=np.array(histStore[ky]['HISTSTORE'][khy]['bins'])
        counts=np.array(histStore[ky]['HISTSTORE'][khy]['counts'])
        xs=xsData[ky]
        norm=1.0*xs*400*1e3/histStore[ky]["NEVENTS"]
        w=counts*norm
        w2=w*w
#         print(histStore[ky]["NEVENTS"],ky,np.sum(w))
        histsSig['hist'].append(w)
        histsSig['w2'].append(w2)
        histsSig['label'].append(utl.getLabelFromTag(ky))

    histsBkg['w2']=np.array(histsBkg['w2'])  
    histsSig['w2']=np.array(histsSig['w2'])  
    histsBkg['bins']=np.array(bins)  
    histsSig['bins']=np.array(bins)  
    odct={}
    odct['bkg']=histsBkg
    odct['sig']=histsSig
    return odct

if plotIdx==2:
    print("Plotting mu pt")
    for khy in [ 'dimuons_mu1_pt', 'dimuons_mu2_pt' ]:
        f,ax=plt.subplots(3,2,figsize=(12,15))
        ax=np.ndarray.flatten(ax)
        for i,cat in enumerate(ALL_CATS):
            CAT=cat
            histStore=data_dict_catWise[CAT]
            odct=getHistograms(histStore,khy,bkgSamples,signalSamples)
            histsSig=odct['sig']
            histsBkg=odct['bkg']
            BINS=histsBkg['bins']
            hep.histplot( histsBkg['hist'], bins=BINS,stack=True, w2=histsBkg['w2'],
                          histtype="fill",label=histsBkg['label'],ax=ax[i]
                        )
            # ax2=ax.twinx()
            hep.histplot( histsSig['hist'], bins=BINS,stack=False, w2=histsSig['w2'],
                          histtype="fill",label=histsSig['label'],ax=ax[i]
                        )
        #     break
        for i,_ in enumerate(ax):
            ax[i].semilogy()
            ax[i].set_title(ALL_CATS[i],fontsize=10)
            ax[i].legend(loc=1,ncol=2,fontsize=10)
        
            ax[i].set_xlim([0,20])
            ax[i].set_ylim([1e1,1e14])
        
            ax[i].tick_params(axis='both',which='both',labelsize=15)
            # ax2.legend(loc=2,ncol=2,fontsize=10)
        #     hep.cms.label("",year="Phase 2",com=14)
            xlbl=khy.replace('dimuons_','').replace("_"," ").replace("mu","$\mu$")
            ax[i].set_xlabel(xlbl,fontsize=12)
        ofname=f"{prefix}/{khy}.png"
        print(ofname)
        f.savefig(ofname,bbox_inches='tight')



if plotIdx==1:
    print("Plotting mu eta")
    for khy in [ 'dimuons_mu1_eta', 'dimuons_mu2_eta' ]:
        f,ax=plt.subplots(3,2,figsize=(12,15))
        ax=np.ndarray.flatten(ax)
        for i,cat in enumerate(ALL_CATS):
            CAT=cat
            histStore=data_dict_catWise[CAT]
            odct=getHistograms(histStore,khy,bkgSamples,signalSamples)
            histsSig=odct['sig']
            histsBkg=odct['bkg']
            BINS=histsBkg['bins']
            hep.histplot( histsBkg['hist'], bins=BINS,stack=True, w2=histsBkg['w2'],
                          histtype="fill",label=histsBkg['label'],ax=ax[i]
                        )
            # ax2=ax.twinx()
            hep.histplot( histsSig['hist'], bins=BINS,stack=False, w2=histsSig['w2'],
                          histtype="fill",label=histsSig['label'],ax=ax[i]
                        )
        #     break
        for i,_ in enumerate(ax):
            ax[i].semilogy()
            ax[i].set_title(ALL_CATS[i],fontsize=10)
            ax[i].legend(loc=1,ncol=2,fontsize=10)
        
            ax[i].set_xlim([-3,3])
            ax[i].set_ylim([1e1,1e14])
        
            ax[i].tick_params(axis='both',which='both',labelsize=15)
            # ax2.legend(loc=2,ncol=2,fontsize=10)
        #     hep.cms.label("",year="Phase 2",com=14)
            xlbl=khy.replace('dimuons_','').replace("_"," ").replace("mu","$\mu$")
            print(xlbl)
            ax[i].set_xlabel(xlbl,fontsize=12)
        ofname=f"{prefix}/{khy}.png"
        print(ofname)
        f.savefig(ofname,bbox_inches='tight')


if plotIdx==3:
    khy='dimuons_mu1_pfRelIso'
    V=0.2
    for khy in [
                'dimuons_mu1_tkRelIso','dimuons_mu2_tkRelIso',
                'dimuons_mu1_pfRelIso','dimuons_mu2_pfRelIso',
                'dimuons_mu1_pfNeutralRelIso','dimuons_mu2_pfNeutralRelIso',
                'dimuons_mu1_pfChargedRelIso','dimuons_mu2_pfChargedRelIso',
                'dimuons_mu1_puppiRelIso','dimuons_mu2_puppiRelIso']:
        f,ax=plt.subplots(3,2,figsize=(15,15))
        ax=np.ndarray.flatten(ax)
        for i,cat in enumerate(ALL_CATS):
            CAT=cat
            histStore=data_dict_catWise[CAT]
            odct=getHistograms(histStore,khy,bkgSamples,signalSamples)
            histsSig=odct['sig']
            histsBkg=odct['bkg']
            BINS=histsBkg['bins']
            maskR= BINS[:-1] > V
            maskL= BINS[:-1] <=V
        
            for j,lbl in enumerate(histsBkg['label']):
                frac=np.sum(histsBkg['hist'][j][maskL]) / np.sum(histsBkg['hist'][j])
                histsBkg['label'][j]=lbl+f" {frac*100:.1f} %"
            for j,lbl in enumerate(histsSig['label']):
                frac=np.sum(histsSig['hist'][j][maskL]) / np.sum(histsSig['hist'][j])
                histsSig['label'][j]=lbl+f" {frac*100:.1f} %"
        
            hep.histplot( histsBkg['hist'], bins=BINS,stack=False, density=True,#w2=histsBkg['w2'],
                          histtype="step",label=histsBkg['label'],ax=ax[i]
                        )
            # ax2=ax.twinx()
            hep.histplot( histsSig['hist'], bins=BINS,stack=False,density=True,# w2=histsSig['w2'],
                          histtype="step",label=histsSig['label'],ax=ax[i]
                        )
        
        #     break
        for i,_ in enumerate(ax):
        #     ax[i].semilogy()
            ax[i].set_title(ALL_CATS[i],fontsize=10)
            ax[i].legend(loc=1,ncol=2,fontsize=13)
        
            ax[i].set_xlim([0.0,3])
        #     ax[i].set_ylim([1e1,1e14])
        
            ax[i].tick_params(axis='both',which='both',labelsize=15)
            xlbl=khy.replace('dimuons_','').replace("_"," ").replace("mu","$\mu$")
            ax[i].set_xlabel(xlbl,fontsize=12)
            ax[i].axvline(V,linestyle='--',color='k',alpha=0.5)
        ofname=f"{prefix}/{khy}.png"
        print(ofname)
        f.savefig(ofname,bbox_inches='tight')

if plotIdx==4:
    V=1.0
    for khy in [
                'dimuons_mu1_tkIsoSum','dimuons_mu2_tkIsoSum',
                'dimuons_mu1_pfIsoSum','dimuons_mu2_pfIsoSum',
                'dimuons_mu1_pfNeutralIsoSum','dimuons_mu2_pfNeutralIsoSum',
                'dimuons_mu1_pfChargedIsoSum','dimuons_mu2_pfChargedIsoSum',
                'dimuons_mu1_puppiIsoSum','dimuons_mu2_puppiIsoSum']:
        f,ax=plt.subplots(3,2,figsize=(15,15))
        ax=np.ndarray.flatten(ax)
        for i,cat in enumerate(ALL_CATS):
            CAT=cat
            histStore=data_dict_catWise[CAT]
            odct=getHistograms(histStore,khy,bkgSamples,signalSamples)
            histsSig=odct['sig']
            histsBkg=odct['bkg']
            BINS=histsBkg['bins']
            maskR= BINS[:-1] > V
            maskL= BINS[:-1] <=V
        
            for j,lbl in enumerate(histsBkg['label']):
                frac=np.sum(histsBkg['hist'][j][maskL]) / np.sum(histsBkg['hist'][j])
                histsBkg['label'][j]=lbl+f" {frac*100:.1f} %"
            for j,lbl in enumerate(histsSig['label']):
                frac=np.sum(histsSig['hist'][j][maskL]) / np.sum(histsSig['hist'][j])
                histsSig['label'][j]=lbl+f" {frac*100:.1f} %"
        
            hep.histplot( histsBkg['hist'], bins=BINS,stack=False, density=True,#w2=histsBkg['w2'],
                          histtype="step",label=histsBkg['label'],ax=ax[i]
                        )
            # ax2=ax.twinx()
            hep.histplot( histsSig['hist'], bins=BINS,stack=False,density=True,# w2=histsSig['w2'],
                          histtype="step",label=histsSig['label'],ax=ax[i]
                        )
        
        #     break
        for i,_ in enumerate(ax):
        #     ax[i].semilogy()
            ax[i].set_title(ALL_CATS[i],fontsize=10)
            ax[i].legend(loc=1,ncol=2,fontsize=13)
        
            ax[i].set_xlim([0.0,20.0])
        
            ax[i].tick_params(axis='both',which='both',labelsize=15)
            xlbl=khy.replace('dimuons_','').replace("_"," ").replace("mu","$\mu$")
            ax[i].set_xlabel(xlbl,fontsize=12)
            ax[i].axvline(V,linestyle='--',color='k',alpha=0.5)
        ofname=f"{prefix}/{khy}.png"
        print(ofname)
        f.savefig(ofname,bbox_inches='tight')


if plotIdx==5:
    f,ax=plt.subplots(3,2,figsize=(15,15))
    ax=np.ndarray.flatten(ax)
    khy='dimuonsDrMuMu'
    V=1.0
    for i,cat in enumerate(ALL_CATS):
        CAT=cat
        histStore=data_dict_catWise[CAT]
        odct=getHistograms(histStore,khy,bkgSamples,signalSamples)
        histsSig=odct['sig']
        histsBkg=odct['bkg']
        BINS=histsBkg['bins']
        maskR= BINS[:-1] > V
        maskL= BINS[:-1] <= V
        
        for j,lbl in enumerate(histsBkg['label']):
            frac=np.sum(histsBkg['hist'][j][maskL]) / np.sum(histsBkg['hist'][j])
            histsBkg['label'][j]=lbl+f" {frac*100:.1f} %"
        for j,lbl in enumerate(histsSig['label']):
            frac=np.sum(histsSig['hist'][j][maskL]) / np.sum(histsSig['hist'][j])
            histsSig['label'][j]=lbl+f" {frac*100:.1f} %"
            
        hep.histplot( histsBkg['hist'], bins=BINS,stack=False, density=True,#w2=histsBkg['w2'],
                      histtype="fill",label=histsBkg['label'],ax=ax[i]
                    )
        # ax2=ax.twinx()
        hep.histplot( histsSig['hist'], bins=BINS,stack=False,density=True,# w2=histsSig['w2'],
                      histtype="step",label=histsSig['label'],ax=ax[i]
                    )
    
    #     break
    for i,_ in enumerate(ax):
    #     ax[i].semilogy()
        ax[i].set_title(ALL_CATS[i],fontsize=10)
        ax[i].legend(loc=1,ncol=2,fontsize=13)
    
        ax[i].set_xlim([0.0,4])
        
        ax[i].tick_params(axis='both',which='both',labelsize=15)
        xlbl=khy.replace('dimuons_','').replace("_"," ").replace("mu","$\mu$")
        ax[i].set_xlabel(xlbl,fontsize=12)
        ax[i].axvline(V,linestyle='--',color='k',alpha=0.5)
    ofname=f"{prefix}/{khy}.png"
    print(ofname)
    f.savefig(ofname,bbox_inches='tight')



if plotIdx==6:
    print("Plotting mu-mu mass")
    for khy in [ 'dimuonMass']:
        f,ax=plt.subplots(3,2,figsize=(12,15))
        ax=np.ndarray.flatten(ax)
        for i,cat in enumerate(ALL_CATS):
            CAT=cat
            histStore=data_dict_catWise[CAT]
            odct=getHistograms(histStore,khy,bkgSamples,signalSamples)
            histsSig=odct['sig']
            histsBkg=odct['bkg']
            BINS=histsBkg['bins']
            hep.histplot( histsBkg['hist'], bins=BINS,stack=True, w2=histsBkg['w2'],
                          histtype="fill",label=histsBkg['label'],ax=ax[i]
                        )
            # ax2=ax.twinx()
            hep.histplot( histsSig['hist'], bins=BINS,stack=False, w2=histsSig['w2'],
                          histtype="fill",label=histsSig['label'],ax=ax[i]
                        )
        #     break
        for i,_ in enumerate(ax):
            ax[i].semilogy()
            ax[i].set_title(ALL_CATS[i],fontsize=10)
            ax[i].legend(loc=1,ncol=2,fontsize=10)
        
            ax[i].set_xlim([0,12])
            ax[i].set_ylim([1e1,1e14])
        
            ax[i].tick_params(axis='both',which='both',labelsize=15)
            # ax2.legend(loc=2,ncol=2,fontsize=10)
        #     hep.cms.label("",year="Phase 2",com=14)
            xlbl=khy.replace('dimuons_','').replace("_"," ").replace("mu","$\mu$")
            ax[i].set_xlabel(xlbl,fontsize=12)
        ofname=f"{prefix}/{khy}.png"
        print(ofname)
        f.savefig(ofname,bbox_inches='tight')





