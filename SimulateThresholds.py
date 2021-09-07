import SyndromeClass as SC
import math
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
import sys


eta,pEst,dx,shapeAspectRatio,timeAspectRatio = int(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),float(sys.argv[5])

numSamples=10
pList = np.linspace(pEst-.009,pEst+.011,15)
dz = dx*shapeAspectRatio
numMeasurements = int(dx*timeAspectRatio)
try:
    saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+".pk",'rb')
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
        S = SC.Syndrome(dz,dx,numMeasurements,p/16,eta)
        S.GenerateErrors()
        S.MatchErrors()
        S.AttemptCorrectionOfErrors()
        ZLogicalError,Z1LogicalError,Z2LogicalError,XLogicallError,TimeLogicalError=S.FindLogicalErrors()
        if TimeLogicalError or XLogicallError or ZLogicalError or Z1LogicalError or Z2LogicalError:
            logicalErrorCounts[count]+=1
        totalCounts[count]+=1
    count+=1
    end=time.time()
    print(end-start)
saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+".pk",'wb')
pickle.dump((logicalErrorCounts,totalCounts,pList),saveFile)
saveFile.close()

