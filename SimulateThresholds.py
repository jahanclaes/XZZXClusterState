import SyndromeClass as SC
import math
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
import sys


eta,pEst,d = int(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3])

numSamples=100
pList = np.linspace(0.0001,2*pEst,15)
try:
    saveFile=open("simulationData3D_"+str(eta)+"_"+str(d)+".pk",'rb')
    logicalErrorCounts,totalCounts,pList=pickle.load(saveFile)
    saveFile.close()
    print("Loaded")
except:
    logicalErrorCounts = [0 for p in pList]
    totalCounts = [0 for p in pList]
count=0
for p in pList:
    start=time.time()
    for i in range(numSamples):
        print(p,i)
        if eta==10000:
            S = SC.Syndrome(3*d,3*d,1,p,eta)
        elif eta==100:
            S = SC.Syndrome(3*d,3*d,d,p,eta)
        else:
            S = SC.Syndrome(d,d,d,p,eta)
        S.GenerateErrors()
        S.MatchErrors()
        S.AttemptCorrectionOfErrors()
        ZLogicalError,XLogicallError=S.FindLogicalErrors()
        if XLogicallError or ZLogicalError:
            logicalErrorCounts[count]+=1
        totalCounts[count]+=1
    count+=1
    end=time.time()
    print(end-start)
saveFile=open("simulationData3D_"+str(eta)+"_"+str(d)+".pk",'wb')
pickle.dump((logicalErrorCounts,totalCounts,pList),saveFile)
saveFile.close()

