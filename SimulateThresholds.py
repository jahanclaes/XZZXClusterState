import SyndromeClass as SC
import math
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
import sys


eta,pEst,dx,clusterType = int(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]),sys.argv[4]

numSamples=10
pList = np.linspace(.0001,2*pEst,15)
if clusterType==XZZX and eta>50:
    dz = dx*3
    numMeasurements = dx*3
else:
    dz = dx
    numMeasurements = dx
try:
    saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+"_"+clusterType+".pk",'rb')
    logicalErrorCounts,totalCounts,pList=pickle.load(saveFile)
    saveFile.close()
    print("Loaded")
except:
    logicalErrorCounts = [[0,0] for p in pList]
    totalCounts = [0 for p in pList]
count=0
for p in pList:
    start=time.time()
    S = SC.Syndrome(dz,dx,numMeasurements,p,eta,clusterType=clusterType)
    for i in range(numSamples):
        print(p,i)
        S.Clear()
        S.GenerateErrors()
        S.MatchErrors()
        S.AttemptCorrectionOfErrors()
        ZLogicalError,XLogicallError=S.FindLogicalErrors()
        if XLogicallError:
            logicalErrorCounts[count][0]+=1
        if ZLogicalError:
            logicalErrorCounts[count][1]+=1
        totalCounts[count]+=1
    count+=1
    end=time.time()
    print(end-start)
saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+"_"+clusterType+".pk",'wb')
pickle.dump((logicalErrorCounts,totalCounts,pList),saveFile)
saveFile.close()

