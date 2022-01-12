import matplotlib
matplotlib.use('TkAgg')
from matplotlib.ticker import NullFormatter,FormatStrFormatter
import matplotlib.pyplot as plt
import math
import numpy as np
import pickle
import sys
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 10

fig =plt.figure(figsize=(11,8))
axes = fig.subplots(2,3)
thresholds = [[0.04613624164653783,  0.06923109529151684,0.0840381810262429],[0.046730046204960074, 0.029242121189152932,0.029167126104762263]]

y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = (.2,.4,.6,.8))

colorList = ['#E54D5C','#208BBA','#94C595','#F7C59F']
f1,f2=.8,1.2

index1=0
for eta in [1,100,10000]:
    index2=0
    for clusterType in ['XZZX','RHG']:
        threshold = thresholds[index2][index1]
        ax = axes[index2,index1]
        if clusterType=="XZZX" and eta>=100:
            dzList = range(5,9)
        else:
            dzList = range(12,100)
        index3=0
        for dz in dzList:
            try:
                if clusterType=="XZZX":
                    saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+".pk",'rb')
                elif clusterType=="RHG":
                    saveFile=open("simulationData3D_RHG_"+str(eta)+"_"+str(dz)+".pk",'rb')
                logicalErrorCounts,totalCounts,pList=pickle.load(saveFile)
                saveFile.close()
                pList = [p*(1+2/eta) for p in pList]
                errorProbList = [logicalErrorCounts[i]/totalCounts[i] for i in range(len(pList))]
                uncertaintyList = [math.sqrt(abs(errorProbList[i]-errorProbList[i]**2))/math.sqrt(totalCounts[i]) for i in range(len(pList)) if pList[i]>f1*threshold and pList[i]<f2*threshold]
                errorProbList = [errorProbList[i] for i in range(len(pList)) if pList[i]>f1*threshold and pList[i]<f2*threshold]
                pList = [pList[i] for i in range(len(pList)) if pList[i]>f1*threshold and pList[i]<f2*threshold]
                if totalCounts[0]>100:
                    if clusterType=='XZZX' and eta>=100:
                        ax.errorbar(pList,errorProbList,uncertaintyList,label="$d_z="+str(3*dz)+"$",color=colorList[index3])
                    else:
                        ax.errorbar(pList,errorProbList,uncertaintyList,label="$d_z="+str(dz)+"$",color=colorList[index3])
                print(eta,clusterType,dz,totalCounts[0])
                index3+=1
            except:
                pass
        ax.text(threshold*.993,.35,"{:0.3f}".format(threshold),horizontalalignment='right',verticalalignment='top')
        ax.plot([threshold,threshold],[0,1],ls='dotted',c='grey')
        ax.set_title(clusterType+" $\\eta="+str(eta)+"$")
        ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,-2))
        ax.set_yscale("log")
        ax.yaxis.set_minor_locator(y_minor)
        ax.set_ylim(.05,.4)
        ax.set_xlim(threshold*f1,threshold*f2)
        if index1!=0:
            ax.yaxis.set_major_formatter(NullFormatter())
            ax.yaxis.set_minor_formatter(NullFormatter())
        else: 
            ax.yaxis.set_minor_formatter(FormatStrFormatter("%0.2f"))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%0.2f"))
            ax.set_ylabel("$p^L$")
        ax.legend(prop={"size":10},loc='lower right')
 
        index2+=1
    ax.tick_params(axis='x', labelsize=10 )
    ax.tick_params(axis='y', labelsize=10 )
    ax.set_xlabel("$p$")
    index1+=1
plt.savefig("thresholdsPhenom.pdf",bbox_inches='tight')
plt.show()
