import sys,os
import json,argparse
import  tabulate
import glob
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

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
        yld=yld[:-1]+f"+/- {nExpected_err:.2f} "
    return yld



finanme='analysis_summary_v3.json'
prefix='results/plots/v3/'

with open(finanme) as f:
    data=json.load(f)

data_perCAT=data['perCAT']
data_perProc=data['perPROC']

#ALL_CATS=list(data_perCAT.keys())
ALL_CATS=['highptMuMuCentral', 'highptMuMuForward', 'MuMuCentral', 'MuMuForeward','rest']
ALL_PROCS=list(data_perProc.keys())

proc='DARK_PHOTON_M1'
print("Summary for Process : ",proc)
tbl=[]
for cat in data_perProc[proc]:
    row=[cat]
    summary= data_perProc[proc][cat]
    nExpected = summary['nExpected']
    nExpected_err = summary['nExpected_err']
    yld=getYieldInHumanReadableString(nExpected,nExpected_err)
    row.append(yld)
    eff=summary["efficiency"]
    row.append(f"{eff:.3e}")
    tbl.append(row)
print(
       tabulate.tabulate(tbl,
                    headers=["Category","Yield","Efficiency"],
                    tablefmt='simple_grid'
                   )
     )


for prc in ALL_PROCS:
    print(prc)
    #if 'DARK' in prc: continue
    if 'QCD' in prc: continue
    cats = ALL_CATS
    sizes = [data_perProc[prc][cat]['nExpected'] for cat in cats]
    E=data_perProc[prc]['inclusive']['nExpected']
    # "{cat} [{data_perProc[prc][cat]['nExpected']*100/(E+1e-9):.1f} %]"
    lbl=[f"{cat} [{data_perProc[prc][cat]['nExpected']*100/(E+1e-9):.1f} %]" for cat in cats]

    f, ax = plt.subplots(figsize=(5,5))
    _=ax.pie(sizes, labels=lbl,#autopct='%1.1f%%',
             textprops={'fontsize': 14})
    # plt.text(-0.5,0.0,prc)
    plt.title(' '.join(prc.split('_')),fontsize=15)

    ofname=f"{prefix}/pie_{prc}.png"
    print(ofname)
    f.savefig(ofname,bbox_inches='tight')

tbl=[]
f,ax=plt.subplots(3,2,figsize=(10,14))
ax=np.ndarray.flatten(ax)
for i,cat in enumerate(data_perCAT):
    x=[]
    y=[]
    x_=[]
    for j,proc in enumerate(data_perCAT[cat]):
        if 'QCD' in proc: continue
        row=[cat]
        summary= data_perCAT[cat][proc]
        nExpected = summary['nExpected']
        nExpected_err = summary['nExpected_err']
        x.append(proc)
        x_.append(f"{j}")
        y.append(nExpected)
    ax[i].barh(x_,y)
    ax[i].set_xlim([1e2,1e16])
    for j,lbl in enumerate(x):
        ax[i].text(1e10,j,lbl,fontsize=10)
#     ax[i].set_yt("")
    ax[i].semilogx()
    ax[i].set_title(cat,fontsize=10)
    ax[i].tick_params(axis='both',which='both',labelsize=12)
    ax[i].set_yticks([])
ofname=f"{prefix}/pie_{prc}.png"
print(ofname)
f.savefig(ofname,bbox_inches='tight')

