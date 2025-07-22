"""


"""


import glob , os
import argparse
import multiprocessing as mp
import random as rd
import time
import subprocess
from threading import Lock

tag='v35p1'
files=glob.glob("may14/*.root")

NUM_CPUS=4
jobsRunning=False
N=0
taskID=0
completedIdx=0
globalFailedCMDS=[]
lock = Lock()
def execute_timeloop():
    timeElapsed=0
    print("Starting the time job !! ")
    global jobsRunning
    while jobsRunning:
        time.sleep(1)
        timeElapsed+=10
        print(f"Time elapsed  : {timeElapsed} secs, running jobs {taskID} / {N} , completed jobs :{completedIdx}",flush=True)
        if timeElapsed >60:
            break
    print("Ending the time job !! ")
def execute_cmd(cmd):
    #value = rd.random()
    #time.sleep(value*3)
    #print("    >>> PRGRESSING  <<<",flush=True)
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()    
    if err.strip()!='':
        print("ERROR FOR CMD : ",cmd)
    rc = proc.returncode
    if rc !=0:
        print(err,flush=True)
        out= cmd + "\nERROR !! \n"+" - x - x"*10+"\n"+err+"\n\n"+" - x - x"*10+"\n"
    #os.system(out)
    return cmd
def progress_callback(out):
    #with lock:
    if True:
        global taskID,completedIdx,globalFailedCMDS,jobsRunning
        completedIdx+=1
        taskID-=1
        if "ERROR" in out:
            globalFailedCMDS.append( completedIdxcompletedIdx.split("ERROR")[0]  )
        if taskID==0:
            jobsRunning=False
    print(f"Fininshed {completedIdx}/{N}. {taskID} remining!\n > {out}" , flush=True)
    

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inputFile" , help="Input files" , default=None)
    parser.add_argument("-n","--ncpu"      , help="N-CPU" , default=2, type=int)
    #parser.add_argument("-p","--printOnly", help="Do not submit the jobs",action='store_true')
    parser.add_argument("--exec", help="execute the cmd",action='store_true')
    parser.add_argument("--execlocal", help="execute the cmd in a sequential queue",action='store_true')
    parser.add_argument("-t","--isTest", help="just process only the one input file",action='store_true')
    args = parser.parse_args()
    
    NUM_CPUS=args.ncpu
    allCm=[]
    if args.inputFile:
        print("Opening input file ",args.inputFile)
        allCmds=[]
        with open(args.inputFile) as f:
            txt=f.readlines()
            for l in txt:
                if l[0]=='#' : continue
                allCmds.append(l[:-1])
    else:
        print("Please give an input file with list of cmds !")

    pool=None
    if args.exec:
        pool = mp.Pool(NUM_CPUS)
    
    proctSet={}
    for i,cmd in enumerate(allCmds):
        print("Procesing ",i," / ",len(allCmds) )
        if args.execlocal:
            print(f"  Executing in sequence !! ")
            os.system(cmd)
        elif args.exec:
            taskID+=1
            N+=1
            print(f"  Adding task to the pool ! n-Task {taskID} \n\t {cmd}")
            pool.apply_async(execute_cmd, args=(cmd,),callback=progress_callback)    
            if taskID==1:
                jobsRunning=True
                pool.apply_async(execute_timeloop)    
        else:
            print(cmd)
        print("= "*20)
        if args.isTest:
            break
    
    
    if args.exec:    
        pool.close()
        pool.join()
        print("Failed cmds")
        for ky in globalFailedCMDS:
            print(ky)
