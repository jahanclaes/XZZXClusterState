import matplotlib
matplotlib.use('TkAgg')
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import pickle
import math
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14


index1 = 0
pEstArray = [[.047, 0.04479042290342928, 0.05076702222377273, 0.06525719184296593, 0.07, 0.07787583933578261, 0.0794256924867034, 0.08287947064179474, 0.08237092892073404],[.047, .034, 0.03188780393893249, 0.03061935675534757, 0.029360441341084065, 0.02913821806346027, 0.029775133390086143, 0.029562543204491906, 0.029518011148030522]]
thresholdsXZZX = []
thresholdsRHG = []
errorsXZZX = []
errorsRHG = []
etaList = [1,5,10,50,100,500,1000,5000,10000]
for eta in etaList:
    index2 = 0
    for clusterType in ['XZZX','RHG']:
        pEst = pEstArray[index2][index1]
        overallPList,overallDList,overallErrorList,overallErrorUncertaintyList = [],[],[],[]
        if eta>=100 and clusterType=="XZZX":
            dzList = range(5,9)
        else:
            dzList = range(12,100)

        for dz in dzList:
            try:
                if clusterType=="XZZX":
                    saveFile=open("simulationData3D_"+str(eta)+"_"+str(dz)+".pk",'rb')
                else:
                    saveFile=open("simulationData3D_RHG_"+str(eta)+"_"+str(dz)+".pk",'rb')
                logicalErrorCounts,totalCounts,pList=pickle.load(saveFile)
                saveFile.close()
                pList = [p*(1+2/eta) for p in pList]
                errorProbList = [logicalErrorCounts[i]/totalCounts[i] for i in range(len(pList))]
                for i in range(len(pList)):
                    if abs(pList[i]-pEst)/pEst<.2:
                        overallPList.append(pList[i])
                        overallDList.append(dz)
                        overallErrorList.append(errorProbList[i])
                        overallErrorUncertaintyList.append(math.sqrt(abs(errorProbList[i]-errorProbList[i]**2))/math.sqrt(totalCounts[i]))
            except:
                pass
        def f(x,pc,nu,A,B,C):
            p,L=x[0],x[1]
            xp = (p-pc)*L**(1/nu)
            return A+B*xp+C*xp**2
        
        overallDataList = np.array([[overallPList[i],overallDList[i]] for i in range(len(overallPList))])
        try:
            popt,pcov=curve_fit(f,overallDataList.T,overallErrorList,p0=[pEst,1,1,1,.1],sigma=overallErrorUncertaintyList)
            print(clusterType,popt[1:])
            pc,nu = popt[0],popt[1]
            pcError,nuError = math.sqrt(pcov[0,0]),math.sqrt(pcov[1,1])
        except:
            pc,nu=float("nan"),float("nan")
            pcError,nuError = float("nan"),float("nan")
        if clusterType=="XZZX":
            thresholdsXZZX.append(pc)
            errorsXZZX.append(pcError)
        if clusterType=="RHG":
            thresholdsRHG.append(pc)
            errorsRHG.append(pcError)
        index2+=1
    index1+=1

print(thresholdsXZZX)
print(thresholdsRHG)
fig = plt.figure(figsize=(4,4))
ax = fig.gca()
ax.errorbar(etaList,thresholdsRHG,errorsRHG,fmt='s',label="RHG",color='#E54D5C')
ax.errorbar(etaList,thresholdsXZZX,errorsXZZX,fmt='^',label="XZZX",color='#208BBA')
ax.set_ylim(0,.1)
ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,-2))
ax.set_xlabel("Bias $\\eta$")
ax.set_ylabel("Threshold error rate $p_{\mathrm{th}}$")
ax.set_xscale("log")
handles, labels = ax.get_legend_handles_labels()
order = [1,0]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
ax.tick_params(axis='x', labelsize=10 )
ax.tick_params(axis='y', labelsize=10 )
plt.savefig("thresholdGraphPhenom.pdf",bbox_inches='tight')
plt.show()
