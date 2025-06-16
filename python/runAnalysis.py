import sys,os
import json,argparse
import  tabulate 
import glob

import util as utl


parser = argparse.ArgumentParser()
parser.add_argument("-m","--mode", help="Run mode",default='analysis')
parser.add_argument("-v","--ver", help="Version")
parser.add_argument("-d","--doc", help="Documentation",default="default selections, Dimuon Analysis")
parser.add_argument("--exec", help="execute the commands",action='store_true')
parser.add_argument("-t","--isTest", help="execute the test command",action='store_true')
parser.add_argument("--isTest2", help="execute the test command",action='store_true')
#args = parser.parse_args()

def run_histMaker():
    cmd_tpl='python3 python/analysis.py -i @@FILE   -t @@TAG    -d @@DEST --doc "@@DOC"'
    test_ext=" -e 100"
    driver_file="data/v1_files.json"
    #driver_file="data/v2_with5parTraking_files.json"
    
    with open(driver_file) as f:
        fileMap=json.load(f)
    
    parser.add_argument("-p","--procs", help="Procsses to process",default=None)
    args = parser.parse_args()
    procs=args.procs
    if procs is not None:
        procs=procs.split("/")
    base=f"results/analysis/{args.ver}/"
    doc=args.doc
        
    for proc in fileMap:
        if procs is not None:
            if proc not in procs: continue 

        for i,fl in enumerate(fileMap[proc]):
    
            destination=base+f'{proc}/parts/'
            tag=proc+f'_{i}'
            if len(fileMap[proc])==1:
                destination=base+f'{proc}/'
                tag=proc
            cmd=str(cmd_tpl)
            cmd=cmd.replace("@@FILE",fl)
            cmd=cmd.replace("@@TAG",tag)
            cmd=cmd.replace("@@DEST",destination)
            cmd=cmd.replace("@@DOC",doc)
            if args.isTest:
               cmd+=test_ext 
            print(cmd)
            if args.exec:
                os.system(cmd)
    
            if args.isTest:
                break
        if args.isTest2:
            break
def run_summaryAnalyzer():
    
    from printSummary import getSummary
    
    parser.add_argument( "--tag", help="Tag",default=None)
    parser.add_argument("-o","--ofname", help="Output summary Json filename",default=None)
    parser.add_argument("-i","--inputFile", help="Input summary Json",default=None)
    args = parser.parse_args()
    if args.inputFile is None:
        base='results/analysis/v2/'
        cats=['highptMuMuCentral','highptMuMuForward','MuMuCentral','MuMuForeward','inclusive','rest']
        cmdExtras=" "

        procFolders=glob.glob(base+'/*')
        procs=[fd.split("/")[-1] for fd in procFolders]
        data_dict_catWise={}
        data_dict_procWise={}
        for cat in cats:
            print("= "*10," Processing category ",cat)
            for proc in procs:
                srchstr=f"{base}/{proc}/*_{cat}_*"
                flist=glob.glob(srchstr)
                if len(flist)!=1:
                    print(f"Ambiguty in the file search ! search string : {srchstr}")
                    exit(1)
                fl=flist[0]
                cmd=f"python3 python/printSummary.py -i {fl} -p {proc}"

                cmd+=cmdExtras
                print(cmd)
                if args.exec:
                    os.system(cmd)
                    with open('summary.json') as f:
                        summary=json.load(f)
                    if cat not in data_dict_catWise:
                        data_dict_catWise[cat]={}
                    if proc not in data_dict_procWise:
                        data_dict_procWise[proc]={}
                    data_dict_catWise[cat][proc] =summary
                    data_dict_procWise[proc][cat]=summary
                    cmd='rm summary.json'
                    os.system(cmd)
            if args.isTest:
                break
            print()
        if args.ofname is not None:
            foutname=args.ofname
        else:
            foutname='analysis_summary.json'
            if args.tag is not None:
                foutname=f'analysis_summary_{args.tag}.json'
        odct={'perCAT':data_dict_catWise,'perPROC':data_dict_procWise}
        with open(foutname,'w') as f:
            json.dump(odct,f,indent=4)
    else:
        print("Loading the summary stats from ",args.inputFile)
        with open(args.inputFile) as f:
            idct=json.load(f)
        data_dict_catWise=idct['perCAT']
        data_dict_procWise=idct['perPROC']
    for cat in data_dict_catWise:
        print("Summary for CAT : ",cat)
        tbl=[]
        for proc in data_dict_catWise[cat]:
            row=[proc]
            summary= data_dict_catWise[cat][proc]
            nExpected = summary['nExpected']
            nExpected_err = summary['nExpected_err']
            yld=utl.getYieldInHumanReadableString(nExpected,nExpected_err)
            row.append(yld)
            eff=summary["efficiency"]
            row.append(f"{eff:.3e}")
            tbl.append(row)
        print(
               tabulate.tabulate(tbl,
                            headers=["Process","Yield","Efficiency"],
                            tablefmt='simple_grid'
                           )
             )

    for proc in data_dict_procWise:
        print("Summary for Process : ",proc)
        tbl=[]
        for cat in data_dict_procWise[proc]:
            row=[cat]
            summary= data_dict_procWise[proc][cat]
            nExpected = summary['nExpected']
            nExpected_err = summary['nExpected_err']
            yld=utl.getYieldInHumanReadableString(nExpected,nExpected_err)
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


def main():
    mode='summary'
    mode='analysis'
    if mode=='analysis':
        run_histMaker()
    if mode=='summary':
        run_summaryAnalyzer()
    #run_summaryAnalyzer()

if __name__=='__main__':
    main()

            
