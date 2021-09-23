import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import numpy as np
import pickle
import sys

#eta,clusterType = int(sys.argv[1]),sys.argv[2]

fig =plt.figure(figsize=(14,11))
axes = fig.subplots(2,3)
index1=0
for eta in [1,100,10000]:
    index2=0
    for clusterType in ['XZZX','RHG']:
        ax = axes[index2,index1]
        dzList = range(5,100)
        if (eta==100 or eta==10000) and clusterType=="XZZX":
            dzList = range(10,100)
        for dz in dzList:
            try:
                saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+"_"+clusterType+".pk",'rb')
                logicalErrorCounts,totalCounts,pList=pickle.load(saveFile)
                print(logicalErrorCounts,totalCounts)
                saveFile.close()
                errorProbList = [sum(logicalErrorCounts[i])/totalCounts[i] for i in range(len(pList))]
                errorXProbList = [logicalErrorCounts[i][0]/totalCounts[i] for i in range(len(pList))]
                errorZProbList = [logicalErrorCounts[i][1]/totalCounts[i] for i in range(len(pList))]
                if totalCounts[0]>500:
                    ax.errorbar(pList,errorProbList,[math.sqrt(abs(errorProbList[i]-errorProbList[i]**2))/math.sqrt(totalCounts[i]) for i in range(len(pList))],label="dz="+str(dz))
            except:
                pass
        ax.set_title(clusterType+" $\\eta="+str(eta)+"$")
        ax.set_xlabel("$p_z$")
        ax.set_ylabel("$p_x^L+p_y^L$")
        ax.set_ylim(0,.5)
        if eta==1:
            ax.set_xlim(.002,.012)
        elif clusterType=="RHG":
            ax.set_xlim(.01,.02)
        elif eta==100 and clusterType=='XZZX': 
            ax.set_xlim(.02,.04)
        elif eta==10000 and clusterType=='XZZX': 
            ax.set_xlim(.03,.06)
        ax.legend()
        index2+=1
    index1+=1
plt.savefig("thresholds.png",bbox_inches='tight')
plt.show()
