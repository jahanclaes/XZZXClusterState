import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import numpy as np
import pickle
import sys

eta = int(sys.argv[1])

if eta==1:
    dzList = list(range(2,9))+[10,12]
if eta==10:
    dzList = [2*i for i in range(2,11)]
if eta==100:
    dzList = [5*i for i in range(2,6)]


for dz in dzList:
    saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+".pk",'rb')
    logicalErrorCounts,totalCounts,pList=pickle.load(saveFile)
    print(logicalErrorCounts,totalCounts)
    saveFile.close()
    errorProbList = [logicalErrorCounts[i]/totalCounts[i] for i in range(len(pList))]
    plt.errorbar(pList,errorProbList,[math.sqrt(errorProbList[i]-errorProbList[i]**2)/math.sqrt(totalCounts[i]) for i in range(len(pList))],label=dz)
plt.legend()
plt.show()
