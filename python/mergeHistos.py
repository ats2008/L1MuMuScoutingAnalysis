import json,glob
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inputFiles", help="inputFile",default='')
    parser.add_argument("-s","--searchString", help=" search string",default=None)
    #parser.add_argument("-t","--tag", help=" tag",default='v0')
    parser.add_argument("-o","--outFile", help="Output json file",default='merjed_output.json')
    args = parser.parse_args()
    
    inputFiles={}
    flist=[]
    if len(args.inputFiles) >0:
        flist+=args.inputFiles.split(",")
    if args.searchString is not None:
        searchString=args.searchString.replace("@","*")
        print("Search string : ",searchString)
        flist += glob.glob(searchString) 
    print("     n-Files  : ",len(flist))
    #inputFiles[args.tag]=flist
    inputFiles['def']=flist
    hist_store={}
    for ky in inputFiles:
        hist_store[ky]=[]
        for fl in inputFiles[ky]:
            if len(fl) < 1: continue
            print("  >  Processing file ",fl)
            with open(fl) as f:
                hist_store[ky].append(json.load(f))
    
    histStore={}
    for ky in hist_store:
        histStore[ky]={"NEVENTS":0,"HISTSTORE":{},"METADATA":{}}
        histStore[ky]["METADATA"]={}
        histStore[ky]["METADATA"]["NMERGE"]=len(hist_store[ky])
        histStore[ky]["METADATA"]["merge_metas"]=[]
        for dta in hist_store[ky]:
            histStore[ky]["METADATA"]["merge_metas"].append(dta["METADATA"])
            histStore[ky]["NEVENTS"]+=dta["NEVENTS"]
            for khy in dta['HISTSTORE']:
                if khy not in histStore[ky]["HISTSTORE"]:
                    histStore[ky]["HISTSTORE"][khy]={}
                    histStore[ky]['HISTSTORE'][khy]['counts']=np.array(dta['HISTSTORE'][khy]['counts'])
                    histStore[ky]['HISTSTORE'][khy]['bins'] =np.array(dta['HISTSTORE'][khy]['bins'])
                    histStore[ky]['HISTSTORE'][khy]['NAME'] =dta['HISTSTORE'][khy]['NAME']
                    histStore[ky]['HISTSTORE'][khy]['DOC'] =dta['HISTSTORE'][khy]['DOC']
                else:
                    histStore[ky]['HISTSTORE'][khy]['counts']+=np.array(dta['HISTSTORE'][khy]['counts'])

                for  c_ in ['counts','bins']:
                        histStore[ky]['HISTSTORE'][khy][c_]=[ float(i) for i in dta['HISTSTORE'][khy][c_] ]

        ofile=args.outFile
        with open(ofile,'w') as f:
            print("Output file ",ofile)
            json.dump(histStore[ky],f,indent=4)
        
        break


if __name__=='__main__':
    main()

            
