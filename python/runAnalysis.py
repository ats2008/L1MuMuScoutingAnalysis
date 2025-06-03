import sys,os
import json,argparse

cmd_tpl='python3 python/analysis.py -i @@FILE   -t @@TAG    -d @@DEST --doc "@@DOC"'
test_ext=" -e 100"
driver_file="data/v1_files.json"

with open(driver_file) as f:
    fileMap=json.load(f)


parser = argparse.ArgumentParser()
parser.add_argument("-v","--ver", help="Version")
parser.add_argument("-d","--doc", help="Documentation",default="default selections, Dimuon Analysis")
#parser.add_argument("--xmin", help=" min bound in histogram",default=-1e6,type=float)
parser.add_argument("--exec", help="execute the commands",action='store_true')
parser.add_argument("--isTest", help="execute the test command",action='store_true')
args = parser.parse_args()

base=f"results/analysis/{args.ver}/"
doc=args.doc
    
for proc in fileMap:
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
