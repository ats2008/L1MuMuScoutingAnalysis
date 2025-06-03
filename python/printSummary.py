import json,glob
import argparse
import numpy as np

XSDBFILE='data/xsdb.json'
xsdb=None
with open(XSDBFILE) as f:
    xsdb=json.load(f)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inputFile", help="inputFile")
    parser.add_argument("-l","--lumi", help=" Target luminosity",default=400.0e3,type=float)
    parser.add_argument("-x","--xs", help=" Input cross-section",default=1.0,type=float)
    parser.add_argument("-n","--hname", help="hitogram to check selections ",default='dimuonMass')
    parser.add_argument("-p","--proc", help="Process name to get xs from",default=None)
    parser.add_argument("--xmin", help=" min bound in histogram",default=-1e6,type=float)
    parser.add_argument("--xmax", help=" max bound in histogram",default=1.0e6,type=float)
    parser.add_argument("-o","--outFile", help="Output fileout",default='merjed_output.json')
    parser.add_argument(  "--printAllHistograms", help="print the details of all the available histograms",action='store_true')
    args = parser.parse_args()
    
    xs=args.xs
    if args.proc is not None:
        if  args.proc not in xsdb:
            print(f"{args.proc} not in xsdb. Available keys are {list(xsdb.keys())}")
            exit()
        else:
            xs=xsdb[args.proc]
    lumi=args.lumi

    hname=args.hname
    xmin=args.xmin
    xmax=args.xmax

    print("  > Input file ",args.inputFile)
    with open( args.inputFile ) as f:
        histStore=json.load(f)
    
    if args.printAllHistograms:
        for i,khy in enumerate(histStore["HISTSTORE"]):
            print(f"  {i+1} ) {khy} : {histStore['HISTSTORE'][khy]['DOC']}" )
        exit(0)
    ntot=histStore['NEVENTS']
    print("  > Number of events            : ",ntot)
    print("  > Selection histogram name    : ",hname)
    dx=histStore["HISTSTORE"][hname]['bins'][1]-histStore["HISTSTORE"][hname]['bins'][0]
    print("  > Selection histogram binning (b0-b1)_edges    : ",dx)
    print("  > Selection histogram selected bounds    : ",f"[{xmin},{xmax}]")
    binx=np.array(histStore["HISTSTORE"][hname]['bins'] )
    mask=np.logical_and(binx< xmax, binx>=xmin)
    counts=np.array(histStore["HISTSTORE"][hname]['counts'])
    npass=np.sum(counts[mask[:-1]])
    print("  > Number of bins in histogram  post selection   : ",np.sum(mask))
    print("  > Number of events in histogram post seelction  : ",npass)
    eff=npass/ntot
    print("        >  selection Efficiency  :   ",f"{eff:.3e}")
    if xs is not None:
        nExpected     = npass/ntot * lumi * xs
        nExpected_err = nExpected/np.sqrt(npass+1e-15)
        if nExpected > 1e8:
            yld=f"{nExpected/1e9:.2f} +/- {nExpected_err/1e9:.2f} B"
        elif nExpected > 1e5:
            yld=f"{nExpected/1e6:.2f} +/- {nExpected_err/1e6:.2f} M"
        elif nExpected > 1e3:
            yld=f"{nExpected/1e6:.2f} +/- {nExpected_err/1e6:.2f} K"
        else :
            yld=f"{nExpected:.2f} +/- {nExpected_err:.2f} K"
            
        print("  target lumi : ",lumi/1e3, " /fb   for xs  = ",xs,"pb Expected  yield : ",yld)
   

if __name__=='__main__':
    main()

            
