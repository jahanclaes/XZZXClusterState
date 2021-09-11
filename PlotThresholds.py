import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import numpy as np
import pickle
import sys

eta,clusterType = int(sys.argv[1]),sys.argv[2]

dzList = range(1,100)
for dz in dzList:
    try:
        saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+"_"+clusterType+".pk",'rb')
        logicalErrorCounts,totalCounts,pList=pickle.load(saveFile)
        print(logicalErrorCounts,totalCounts)
        saveFile.close()
        errorProbList = [sum(logicalErrorCounts[i])/totalCounts[i] for i in range(len(pList))]
        errorXProbList = [logicalErrorCounts[i][0]/totalCounts[i] for i in range(len(pList))]
        errorZProbList = [logicalErrorCounts[i][1]/totalCounts[i] for i in range(len(pList))]
        if totalCounts[0]>200:
            #plt.errorbar(pList,errorProbList,[math.sqrt(errorProbList[i]-errorProbList[i]**2)/math.sqrt(totalCounts[i]) for i in range(len(pList))],label=dz)
            plt.subplot(211)
            plt.title("X")
            plt.errorbar(pList,errorXProbList,[math.sqrt(errorXProbList[i]-errorXProbList[i]**2)/math.sqrt(totalCounts[i]) for i in range(len(pList))],label=dz)
            plt.subplot(212)
            plt.title("Z")
            plt.errorbar(pList,errorZProbList,[math.sqrt(errorZProbList[i]-errorZProbList[i]**2)/math.sqrt(totalCounts[i]) for i in range(len(pList))],label=dz)
    except:
        pass
plt.legend()
plt.show()
