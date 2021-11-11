import matplotlib
#matplotlib.use('TkAgg')
import numpy as np
import networkx as nx
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import random
import ctypes
import pickle
import time

# load pypm library
# declare function argument types
_libpypm = ctypes.cdll.LoadLibrary("./libpypm.so")
_libpypm.infty.argtypes = None
_libpypm.infty.restype = ctypes.c_int
_libpypm.mwpm.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
_libpypm.mwpm.restype = None

INFTY = _libpypm.infty() # Integer that represents infinity (Blossom V algorithm). (Can be useful when converting float to int weights for MWPM)

PauliList = ["I","X","Y","Z"]
PauliArray = [["I","X","Y","Z"],["X","I","Z","Y"],["Y","Z","I","X"],["Z","Y","X","I"]]

def PauliMultiply(P1,P2):
    return PauliArray[PauliList.index(P1)][PauliList.index(P2)]

def Hadamard(P1):
    if P1=="X":
        return "Z"
    elif P1=="Y":
        return "Y"
    elif P1=="Z":
        return "X"
    elif P1=="I":
        return "I"

def mwpm_ids(edges):
    """Minimum Weight Perfect Matching using node ids (Blossom V algorithm).

    Notes:

    * Node ids are assumed to form a contiguous set of non-negative integers starting at zero, e.g.  {0, 1, ...}.
    * All nodes are assumed to participate in at least one edge.

    :param edges: Edges as [(node_id, node_id, weight), ...].
    :type edges: list of (int, int, int)
    :return: Set of matches as {(node_id, node_id), ...}. (Each tuple is sorted.)
    :rtype: set of (int, int)
    """
    # extract and sort node ids
    node_ids = sorted(set(id for (id_a, id_b, _) in edges for id in (id_a, id_b)))
    # count n_nodes
    n_nodes = len(node_ids)
    # check node ids form contiguous set of non-negative integers starting at zero
    assert n_nodes == 0 or (node_ids[0] == 0 and node_ids[-1] == n_nodes - 1), (
        'Node ids are not a contiguous set of non-negative integers starting at zero.')
    # count n_edges
    n_edges = len(edges)
    # unzip edges
    nodes_a, nodes_b, weights = zip(*edges) if n_edges else ([], [], [])
    # prepare array types
    mates_array_type = ctypes.c_int * n_nodes
    edges_array_type = ctypes.c_int * n_edges
    # prepare empty mates
    mates_array = mates_array_type()
    # call C interface
    _libpypm.mwpm(ctypes.c_int(n_nodes), mates_array, ctypes.c_int(n_edges),edges_array_type(*nodes_a), edges_array_type(*nodes_b), edges_array_type(*weights))
    # structure of mates: mates_array[i] = j means ith node matches jth node
    # convert to more useful format: e.g. convert [1, 0, 3, 2] to {(0, 1), (2, 3)}
    mates = {tuple(sorted((a, b))) for a, b in enumerate(mates_array)}
    return mates


def mwpm(edges):
    """Minimum Weight Perfect Matching using node objects (Blossom V algorithm).

    :param edges: Edges as [(node, node, weight), ...].
    :type edges: list of (object, object, int)
    :return: Set of matches as {(node, node), ...}.
    :rtype: set of (object, object)
    """
    # list of nodes without duplicates
    nodes = list(set(node for (node_a, node_b, _) in edges for node in (node_a, node_b)))
    # dict of node to id
    node_to_id = dict((n, i) for i, n in enumerate(nodes))
    # edges using ids
    edge_ids = [(node_to_id[node_a], node_to_id[node_b], weight) for node_a, node_b, weight in edges]
    # mwpm using ids
    mate_ids = mwpm_ids(edge_ids)
    # matches using objects
    mates = {(nodes[node_id_a], nodes[node_id_b]) for node_id_a, node_id_b in mate_ids}
    return mates



"""
Overall measurement pattern on a single stabilizer is X-X-Z-Z:
  Z      2
 / \    / \ 
X   X  1   4 
 \ /    \ /
  Z      3  
"""

class Syndrome:

    def __init__(self,dz,dx,dt,p,eta,clusterType='XZZX'):
        """
        Initializes a Syndrome object
        dx and dz are the dimensions of the code
        dt is the number of rounds of stabilizer measurements we are simulating
        """
        self.dx=dx
        self.dz=dz
        self.dt = dt
        self.clusterType=clusterType
        self.correctMatches = {}
        self.decodedMatches = {}
        self.defectDictionary = {}
        self.numCrossingsInTime = 0
        pIdle=p
        pPrep=p
        pMeasure=p
        pCNOT=p
        pCZ=p
        singleQubitErrors = ["I","Z","X","Y"]
        twoQubitErrors =["II","IZ","ZI","ZZ","IX","IY","XI","XX","XY","XZ","YI","YX","YY","YZ","ZX","ZY"] 
        weightsSPAM = [1-(pMeasure+pPrep)*(1+2/eta),pMeasure+pPrep,(pMeasure+pPrep)/eta,(pMeasure+pPrep)/eta]
        weightsIdle = [0]*4 
        weightsCNOT = [1-2*pCNOT-12*pCNOT/eta,pCNOT,pCNOT/2,pCNOT/2]+[pCNOT/eta]*12 
        weightsCZ = [1-2*pCZ-pCZ**2-12*pCZ/eta,pCZ,pCZ,pCZ**2]+[pCZ/eta]*12
        self.SPAMErrors = {singleQubitErrors[i]:weightsSPAM[i] for i in range(4)}
        self.idleErrors = {singleQubitErrors[i]:weightsIdle[i] for i in range(4)}
        if clusterType=="XZZX":
            self.CNOTErrors = {twoQubitErrors[i]:weightsCNOT[i] for i in range(16)} #Errors are (target,control)
        else:
            self.CNOTErrors = {Hadamard(twoQubitErrors[i][0])+twoQubitErrors[i][1]:weightsCZ[i] for i in range(16)} #Hadamards make this the RHG Lattice
        self.CZErrors = {twoQubitErrors[i]:weightsCZ[i] for i in range(16)}
        try:
            saveFile = open("distanceDict_"+str(clusterType)+"_"+str(dt)+"_"+str(dz)+"_"+str(dx)+"_"+str(eta)+".pk",'rb')
            self.distanceDict,baseProbability = pickle.load(saveFile)
            saveFile.close()
            for key1 in self.distanceDict:
                for key2 in self.distanceDict[key1]:
                    self.distanceDict[key1][key2]=self.distanceDict[key1][key2]-math.log(baseProbability/p)
            print("Loaded")
        except:
            self.distanceDict = self.GenerateDistanceDict()
            saveFile = open("distanceDict_"+str(clusterType)+"_"+str(dt)+"_"+str(dz)+"_"+str(dx)+"_"+str(eta)+".pk",'wb')
            pickle.dump((self.distanceDict,p),saveFile)
            saveFile.close()

   
    def Clear(self):
        self.correctMatches = {}
        self.decodedMatches = {}

    def PlotCluster(self,ax,plotScaffold=True):
        #Plot the cluster qubits
        if plotScaffold:
            clusterColor,clusterLineWeight="grey",1
            ZQubitSize,XQubitSize=100,10
            for tIndex in range(self.dt):
                for zIndex in range(2*self.dz-1):
                    for xIndex in range(self.dx-zIndex%2):
                        if zIndex%2==0:
                            ax.scatter(zIndex/2,xIndex+(zIndex%2)/2,tIndex,marker='o',color='k',s=XQubitSize)
                            ax.scatter(zIndex/2,xIndex+(zIndex%2)/2,tIndex+1/2,marker='$\circ$',color='k',s=ZQubitSize)
                        else:
                            ax.scatter(zIndex/2,xIndex+(zIndex%2)/2,tIndex,marker='$\circ$',color='k',s=ZQubitSize)
                            ax.scatter(zIndex/2,xIndex+(zIndex%2)/2,tIndex+1/2,marker='o',color='k',s=XQubitSize)
                    for xIndex in range(self.dx-(zIndex+1)%2):
                        if zIndex%2==0:
                            ax.scatter(zIndex/2,xIndex+((zIndex+1)%2)/2,tIndex,marker='o',color='k',s=XQubitSize)
                        else:
                            ax.scatter(zIndex/2,xIndex+((zIndex+1)%2)/2,tIndex+1/2,marker='o',color='k',s=XQubitSize)
                for zIndex in range(2*self.dz-1):
                    plt.plot([zIndex/2,zIndex/2],[0,self.dx-1],[tIndex+(zIndex%2)/2,tIndex+(zIndex%2)/2],color=clusterColor,lw=clusterLineWeight)
                for xIndex in range(self.dx):
                    plt.plot([0,self.dz-1],[xIndex,xIndex],[tIndex+1/2,tIndex+1/2],color=clusterColor,lw=clusterLineWeight)
                    if xIndex<self.dx-1:
                        plt.plot([0,self.dz-1],[xIndex+1/2,xIndex+1/2],[tIndex,tIndex],color=clusterColor,lw=clusterLineWeight)
            for zIndex in range(2*self.dz-1):
                for xIndex in range(self.dx-(zIndex)%2):
                    plt.plot([zIndex/2,zIndex/2],[xIndex+(zIndex%2)/2,xIndex+(zIndex%2)/2],[0,self.dt-1/2],color=clusterColor,lw=clusterLineWeight)
        #Plot the errors
        errorColor,errorWeight,matchColor='red',4,'blue'
        for key in self.correctMatches:
            point1,point2 = key, self.correctMatches[key]
            if type(point1)==str:
                point1,point2=point2,point1
            if point2=="T" and type(point1)!=str:
                point2 = (point1[0],point1[1],self.dx-1)
            elif point2=="B" and type(point1)!=str:
                point2 = (point1[0],point1[1],-1)
            elif point2=="L" and type(point1)!=str:
                point2 = (point1[0],-1,point1[2])
            elif point2=="R" and type(point1)!=str:
                point2 = (point1[0],2*self.dz-1,point1[2])
            if type(point1)!=str:
                if point1[1]%2==0:
                    plt.plot([point1[1]/2,point2[1]/2],[point1[2]+((point1[1]+1)%2)/2,point2[2]+((point2[1]+1)%2)/2],[point1[0]+.5,point2[0]+.5],color=errorColor,lw=errorWeight)
                else:
                    plt.plot([point1[1]/2,point2[1]/2],[point1[2]+((point1[1]+1)%2)/2,point2[2]+((point2[1]+1)%2)/2],[point1[0]+1,point2[0]+1],color=errorColor,lw=errorWeight)
        for point1,point2 in self.decodedMatches:
            if type(point1[0])==str:
                point1,point2=point2,point1
            if point2[0]=="T" and type(point1[0])!=str:
                point2 = (point1[0],point1[1],self.dx-1)
            elif point2[0]=="B" and type(point1[0])!=str:
                point2 = (point1[0],point1[1],-1)
            elif point2[0]=="L" and type(point1[0])!=str:
                point2 = (point1[0],-1,point1[2])
            elif point2[0]=="R" and type(point1[0])!=str:
                point2 = (point1[0],2*self.dz-1,point1[2])
            if type(point1[0])!=str:
                if point1[1]%2==0:
                    plt.plot([point1[1]/2,point2[1]/2],[point1[2]+((point1[1]+1)%2)/2,point2[2]+((point2[1]+1)%2)/2],[point1[0]+.5,point2[0]+.5],color=matchColor,lw=errorWeight,ls='dotted')
                else:
                    plt.plot([point1[1]/2,point2[1]/2],[point1[2]+((point1[1]+1)%2)/2,point2[2]+((point2[1]+1)%2)/2],[point1[0]+1,point2[0]+1],color=matchColor,lw=errorWeight,ls='dotted')

    def AddMatchedPair(self,point1,point2):
        if type(point1)!=str and type(point2)!=str and abs(point1[0]-point2[0])>self.dt-abs(point1[0]-point2[0]):
            self.numCrossingsInTime+=1
        if point1!=point2:
            if point1 in self.correctMatches.keys() and point2 in self.correctMatches.keys():
                if self.correctMatches[point1]!=self.correctMatches[point2]:
                    self.correctMatches[self.correctMatches[point1]]=self.correctMatches[point2]
                    self.correctMatches[self.correctMatches[point2]]=self.correctMatches[point1]
                    del self.correctMatches[point1]
                    del self.correctMatches[point2]
                else:
                    del self.correctMatches[self.correctMatches[point1]]
                    del self.correctMatches[point1]
                    del self.correctMatches[point2]
            elif point1 in self.correctMatches.keys():
                self.correctMatches[self.correctMatches[point1]]=point2
                self.correctMatches[point2]=self.correctMatches[point1]
                del self.correctMatches[point1]
            elif point2 in self.correctMatches.keys():
                self.correctMatches[self.correctMatches[point2]]=point1
                self.correctMatches[point1]=self.correctMatches[point2]
                del self.correctMatches[point2]
            else:
                self.correctMatches[point1]=point2
                self.correctMatches[point2]=point1
            for key in self.correctMatches.keys():
                if self.correctMatches[self.correctMatches[key]]!=key:
                    print("Error in reflexivity",key,self.correctMatches[key],self.correctMatches[self.correctMatches[key]])
                if type(key)!=str and type(self.correctMatches[key])!=str and (key[1]+self.correctMatches[key][1])%2!=0:
                    print("Error in parity",key,self.correctMatches[key])
                    print(asdf)

    def AddMatchedPairDefectDict(self,point1,point2):
        if point1!=point2:
            if point1 in self.defectDictionary.keys() and point2 in self.defectDictionary.keys():
                if self.defectDictionary[point1]!=self.defectDictionary[point2]:
                    self.defectDictionary[self.defectDictionary[point1]]=self.defectDictionary[point2]
                    self.defectDictionary[self.defectDictionary[point2]]=self.defectDictionary[point1]
                    del self.defectDictionary[point1]
                    del self.defectDictionary[point2]
                else:
                    del self.defectDictionary[self.defectDictionary[point1]]
                    del self.defectDictionary[point1]
                    del self.defectDictionary[point2]
            elif point1 in self.defectDictionary.keys():
                self.defectDictionary[self.defectDictionary[point1]]=point2
                self.defectDictionary[point2]=self.defectDictionary[point1]
                del self.defectDictionary[point1]
            elif point2 in self.defectDictionary.keys():
                self.defectDictionary[self.defectDictionary[point2]]=point1
                self.defectDictionary[point1]=self.defectDictionary[point2]
                del self.defectDictionary[point2]
            else:
                self.defectDictionary[point1]=point2
                self.defectDictionary[point2]=point1
            for key in self.defectDictionary.keys():
                if self.defectDictionary[self.defectDictionary[key]]!=key:
                    print("Error in reflexivity",key,self.defectDictionary[key],self.defectDictionary[self.defectDictionary[key]])
                if type(key)!=str and type(self.defectDictionary[key])!=str and (key[1]+self.defectDictionary[key][1])%2!=0:
                    print("Error in parity",key,self.defectDictionary[key])
                    print(asdf)


    def GenerateDefectsFromError(self,errorIndex,errorString,tIndex,zIndex,xIndex):
        """
        Generates pairs of defects associated to an error
        ErrorIndex indexes the type of error that we are experiencing
        tIndex,xIndex,zIndex are coordinates of the ancilla measurement qubit, layer-1 data qubits, or layer-2 data qubits
        """
        defectPairList=[] 
        if errorIndex==0: # Ancilla Preparation/Measurement errors
            if errorString == "Z" or errorString=="Y":
                defectPairList.append(((tIndex-1,zIndex,xIndex),(tIndex,zIndex,xIndex)))
        if errorIndex==1: # Layer1 data qubit Preparation/Measurement errors
            if zIndex%2==0 and (errorString == "Z" or errorString=="Y"):
                defectPairList.append(((tIndex-1,zIndex+1,xIndex),(tIndex-1,zIndex-1,xIndex)))
            if zIndex%2==1 and (errorString == "X" or errorString=="Y"):
                defectPairList.append(((tIndex-1,zIndex,xIndex+zIndex%2-1),(tIndex-1,zIndex,xIndex+zIndex%2)))
        if errorIndex==2: # Layer2 data qubit Preparation/Measurement errors
            if zIndex%2==0 and (errorString == "X" or errorString=="Y"):
                defectPairList.append(((tIndex,zIndex,xIndex+zIndex%2-1),(tIndex,zIndex,xIndex+zIndex%2)))
            if zIndex%2==1 and (errorString == "Z" or errorString=="Y"):
                defectPairList.append(((tIndex,zIndex+1,xIndex),(tIndex,zIndex-1,xIndex)))
        if errorIndex==3 and zIndex>0: #Gate to the left of the ancilla
            if errorString[0]=="X" or errorString[0]=="Y":
                defectPairList.append(((tIndex-1+zIndex%2,zIndex-1,xIndex-zIndex%2),(tIndex-1+zIndex%2,zIndex-1,xIndex-zIndex%2+1)))
            if errorString[0]=="Z" or errorString[0]=="Y":
                defectPairList.append(((tIndex-1,zIndex-2,xIndex),(tIndex,zIndex,xIndex)))
            if errorString[1]=="Z" or errorString[1]=="Y":
                defectPairList.append(((tIndex-1,zIndex,xIndex),(tIndex,zIndex,xIndex)))
            if errorString[1]=="X" or errorString[1]=="Y":
                defectPairList.append(((tIndex-1+zIndex%2,zIndex-1,xIndex-zIndex%2),(tIndex-1+zIndex%2,zIndex-1,xIndex-zIndex%2+1)))
        if errorIndex==4 and zIndex<2*self.dz-2: #Gate to the right of the ancilla
            if errorString[0]=="X" or errorString[0]=="Y":
                defectPairList.append(((tIndex-1+zIndex%2,zIndex+1,xIndex-zIndex%2),(tIndex-1+zIndex%2,zIndex+1,xIndex-zIndex%2+1)))
            if errorString[0]=="Z" or errorString[0]=="Y":
                defectPairList.append(((tIndex,zIndex,xIndex),(tIndex,zIndex+2,xIndex)))
            if errorString[1]=="Z" or errorString[1]=="Y":
                defectPairList.append(((tIndex-1,zIndex,xIndex),(tIndex,zIndex,xIndex)))
            if errorString[1]=="X" or errorString[1]=="Y":
                defectPairList.append(((tIndex-1+zIndex%2,zIndex+1,xIndex-zIndex%2),(tIndex-1+zIndex%2,zIndex+1,xIndex-zIndex%2+1)))
                defectPairList.append(((tIndex-1+zIndex%2,zIndex-1,xIndex-zIndex%2),(tIndex-1+zIndex%2,zIndex-1,xIndex-zIndex%2+1)))
        if errorIndex==5 and xIndex>zIndex%2-1: #Gate below the ancilla
            if errorString[0]=="X" or errorString[0]=="Y":
                defectPairList.append(((tIndex-1,zIndex,xIndex-1),(tIndex,zIndex,xIndex)))
            if errorString[0]=="Z" or errorString[0]=="Y":
                defectPairList.append(((tIndex-1+zIndex%2,zIndex-1,xIndex-zIndex%2),(tIndex-1+zIndex%2,zIndex+1,xIndex-zIndex%2)))
            if errorString[1]=="Z" or errorString[1]=="Y":
                defectPairList.append(((tIndex-1,zIndex,xIndex),(tIndex,zIndex,xIndex)))
            if errorString[1]=="X" or errorString[1]=="Y":
                defectPairList.append(((tIndex-1+zIndex%2,zIndex-1,xIndex+1-zIndex%2),(tIndex-1+zIndex%2,zIndex+1,xIndex+1-zIndex%2)))
        if errorIndex==6 and xIndex<self.dx-1: #Gate above the ancilla
            if errorString[0]=="X" or errorString[0]=="Y":
                defectPairList.append(((tIndex,zIndex,xIndex+1),(tIndex,zIndex,xIndex)))
            if errorString[0]=="Z" or errorString[0]=="Y":
                defectPairList.append(((tIndex-1+zIndex%2,zIndex-1,xIndex+1-zIndex%2),(tIndex-1+zIndex%2,zIndex+1,xIndex+1-zIndex%2)))
            if errorString[1]=="Z" or errorString[1]=="Y":
                defectPairList.append(((tIndex-1,zIndex,xIndex),(tIndex,zIndex,xIndex)))
            if errorString[1]=="X" or errorString[1]=="Y":
                pass
        if errorIndex==7: #First Layer Entangling Gate
            if zIndex%2==0:
                if errorString[0]=="X" or errorString[0]=="Y": #Upper end
                    defectPairList.append(((tIndex,zIndex,xIndex+zIndex%2-1),(tIndex,zIndex,xIndex+zIndex%2)))
                if errorString[0]=="Z" or errorString[0]=="Y": #Upper end
                    defectPairList.append(((tIndex-1,zIndex-1,xIndex),(tIndex-1,zIndex+1,xIndex)))
                if errorString[1]=="Z" or errorString[1]=="Y": #Lower end
                    defectPairList.append(((tIndex-1,zIndex-1,xIndex),(tIndex-1,zIndex+1,xIndex)))
                if errorString[1]=="X" or errorString[1]=="Y": #Lower end
                    pass
            if zIndex%2==1:
                if errorString[0]=="X" or errorString[0]=="Y": #Lower end
                    defectPairList.append(((tIndex-1,zIndex,xIndex+zIndex%2-1),(tIndex-1,zIndex,xIndex+zIndex%2)))  
                if errorString[0]=="Z" or errorString[0]=="Y": #Lower end           
                    pass 
                if errorString[1]=="Z" or errorString[1]=="Y": #Upper end
                    defectPairList.append(((tIndex,zIndex-1,xIndex),(tIndex,zIndex+1,xIndex)))  
                if errorString[1]=="X" or errorString[1]=="Y": #Upper end
                    defectPairList.append(((tIndex-1,zIndex,xIndex+zIndex%2-1),(tIndex-1,zIndex,xIndex+zIndex%2)))    
        if errorIndex==8: #Second Layer Entangling Gate
            if zIndex%2==0:
                if errorString[0]=="X" or errorString[0]=="Y": #Lower end
                    defectPairList.append(((tIndex,zIndex,xIndex+zIndex%2-1),(tIndex,zIndex,xIndex+zIndex%2))) 
                if errorString[0]=="Z" or errorString[0]=="Y": #Lower end           
                    pass 
                if errorString[1]=="Z" or errorString[1]=="Y": #Upper end
                    defectPairList.append(((tIndex,zIndex-1,xIndex),(tIndex,zIndex+1,xIndex)))
                if errorString[1]=="X" or errorString[1]=="Y": #Upper end
                    defectPairList.append(((tIndex,zIndex,xIndex+zIndex%2-1),(tIndex,zIndex,xIndex+zIndex%2)))
            if zIndex%2==1:
                if errorString[0]=="X" or errorString[0]=="Y": #Upper end
                    defectPairList.append(((tIndex,zIndex,xIndex+zIndex%2-1),(tIndex,zIndex,xIndex+zIndex%2)))
                if errorString[0]=="Z" or errorString[0]=="Y": #Upper end
                    defectPairList.append(((tIndex,zIndex-1,xIndex),(tIndex,zIndex+1,xIndex)))
                if errorString[1]=="Z" or errorString[1]=="Y": #Lower end
                    defectPairList.append(((tIndex,zIndex-1,xIndex),(tIndex,zIndex+1,xIndex)))
                if errorString[1]=="X" or errorString[1]=="Y": #Lower end
                    pass
        # Relabel boundary points appropriately
        for index in range(len(defectPairList)):
            point1,point2 = defectPairList[index]
            if point1[0]==-1:
                point1==(self.dt-1,point1[1],point1[2])
            if point1[1]==-1:
                point1="L"
            elif point1[1]==2*self.dz-1:
                point1="R"
            elif point1[2]==-1:
                point1="B"
            elif point1[2]==self.dx-1 and point1[1]%2==0:
                point1="T"
            if point2[0]==-1:
                point2==(self.dt-1,point2[1],point2[2])
            if point2[1]==-1:
                point2="L"
            elif point2[1]==2*self.dz-1:
                point2="R"
            elif point2[2]==-1:
                point2="B"
            elif point2[2]==self.dx-1 and point2[1]%2==0:
                point2="T"
            defectPairList[index]=(point1,point2)
        # Transfer to dictionary and back to list, to get rid of excess pairs
        self.defectDictionary = {} #Clear the defectDictionary
        for point1,point2 in defectPairList:
            self.AddMatchedPairDefectDict(point1,point2)
        defectPairList = []
        for key in self.defectDictionary:
            if not (key in [d[0] for d in defectPairList] or key in [d[1] for d in defectPairList]):
                defectPairList.append((key,self.defectDictionary[key]))
        return defectPairList

    def GenerateErrors(self):
        """
        Generate defects according to some error model.
        Populates self.correctMatches
        """
        for tIndex in range(self.dt):
            for zIndex in range(2*self.dz-1):
                defectPairList = []
                for xIndex in range(self.dx-1+zIndex%2): #Index over ancilla qubits 
                    #Ancilla Preparation/Measurement errors
                    ancillaError = random.choices(list(self.SPAMErrors.keys()),weights=[self.SPAMErrors[key] for key in self.SPAMErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(0,ancillaError,tIndex,zIndex,xIndex)
                    #Gate to the left of the Ancilla
                    gateError = random.choices(list(self.CNOTErrors.keys()),weights=[self.CNOTErrors[key] for key in self.CNOTErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(3,gateError,tIndex,zIndex,xIndex)
                    #Gate to the right of the Ancilla
                    gateError = random.choices(list(self.CNOTErrors.keys()),weights=[self.CNOTErrors[key] for key in self.CNOTErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(4,gateError,tIndex,zIndex,xIndex)
                    #Gate below the Ancilla
                    gateError = random.choices(list(self.CZErrors.keys()),weights=[self.CZErrors[key] for key in self.CZErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(5,gateError,tIndex,zIndex,xIndex)
                    #Gate above the Ancilla
                    gateError = random.choices(list(self.CZErrors.keys()),weights=[self.CZErrors[key] for key in self.CZErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(6,gateError,tIndex,zIndex,xIndex)
                for xIndex in range(self.dx-zIndex%2): #Index over data qubits
                    #Layer 1 Data Preparation/Measurement errors
                    layerOneDataError = random.choices(list(self.SPAMErrors.keys()),weights=[self.SPAMErrors[key] for key in self.SPAMErrors.keys()])[0]
                    if self.clusterType=="RHG" and zIndex%2==1:
                        layerOneDataError = Hadamard(layerOneDataError)
                    defectPairList = defectPairList+self.GenerateDefectsFromError(1,layerOneDataError,tIndex,zIndex,xIndex)
                    #Layer 2 Data Preparation/Measurement errors
                    layerTwoDataError = random.choices(list(self.SPAMErrors.keys()),weights=[self.SPAMErrors[key] for key in self.SPAMErrors.keys()])[0]
                    if self.clusterType=="RHG" and zIndex%2==0:
                        layerTwoDataError = Hadamard(layerTwoDataError)
                    defectPairList = defectPairList+self.GenerateDefectsFromError(2,layerTwoDataError,tIndex,zIndex,xIndex)
                    #First Layer Entangling Gate
                    gateError = random.choices(list(self.CNOTErrors.keys()),weights=[self.CNOTErrors[key] for key in self.CNOTErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(7,gateError,tIndex,zIndex,xIndex)
                    #Second Layer Entangling Gate
                    gateError = random.choices(list(self.CNOTErrors.keys()),weights=[self.CNOTErrors[key] for key in self.CNOTErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(8,gateError,tIndex,zIndex,xIndex)
                for defect1,defect2 in defectPairList:
                    self.AddMatchedPair(defect1,defect2)
                        
    def GenerateDistanceDict(self):
        """
        Generates a lookup table to find the distance between defects.
        Uses GenerateDefectsFromError to find the distances
        Currectly assumes independent probabilities just add--this is a first-order approximation, can modify later
        """
        distanceGraph = nx.Graph()
        def AddToNodeGraph(probability,defectPairList):
            if len(defectPairList)==1 and defectPairList[0][0]!=defectPairList[0][1]:
                defect1,defect2 = defectPairList[0][0],defectPairList[0][1]
                try:
                    distanceGraph[defect1][defect2]['weight']=distanceGraph[defect1][defect2]['weight']+probability
                except:
                    distanceGraph.add_edge(defect1,defect2,weight=probability)
        for tIndex in range(self.dt):
            for zIndex in range(2*self.dz-1):
                defectPairList = []
                for xIndex in range(self.dx-1+zIndex%2): #Index over ancilla qubits 
                    #Ancilla Preparation/Measurement errors
                    for ancillaError in self.SPAMErrors:
                        defectPairList = self.GenerateDefectsFromError(0,ancillaError,tIndex,zIndex,xIndex)
                        probability = self.SPAMErrors[ancillaError]
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
                    #Gate to the left of the Ancilla
                    for gateError in self.CNOTErrors:
                        defectPairList = self.GenerateDefectsFromError(3,gateError,tIndex,zIndex,xIndex)
                        probability = self.CNOTErrors[gateError]
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
                    #Gate to the right of the Ancilla
                    for gateError in self.CNOTErrors:
                        defectPairList = self.GenerateDefectsFromError(4,gateError,tIndex,zIndex,xIndex)
                        probability = self.CNOTErrors[gateError]
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
                    #Gate below the Ancilla
                    for gateError in self.CZErrors:
                        defectPairList = self.GenerateDefectsFromError(5,gateError,tIndex,zIndex,xIndex)
                        probability = self.CZErrors[gateError]
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
                    #Gate above the Ancilla
                    for gateError in self.CZErrors:
                        defectPairList = self.GenerateDefectsFromError(6,gateError,tIndex,zIndex,xIndex)
                        probability = self.CZErrors[gateError]
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
                for xIndex in range(self.dx-zIndex%2): #Index over data qubits
                    #Layer 1 Data Preparation/Measurement errors
                    for layerOneDataError in self.SPAMErrors:
                        probability=self.SPAMErrors[layerOneDataError]
                        if self.clusterType=="RHG" and zIndex%2==1:
                            layerOneDataError = Hadamard(layerOneDataError)
                        defectPairList = self.GenerateDefectsFromError(1,layerOneDataError,tIndex,zIndex,xIndex)
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
                    #Layer 2 Data Preparation/Measurement errors
                    for layerTwoDataError in self.SPAMErrors:
                        probability=self.SPAMErrors[layerTwoDataError]
                        if self.clusterType=="RHG" and zIndex%2==0:
                            layerTwoDataError = Hadamard(layerTwoDataError)
                        defectPairList = self.GenerateDefectsFromError(2,layerTwoDataError,tIndex,zIndex,xIndex)
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
                    #First Layer Entangling Gate
                    for gateError in self.CNOTErrors:
                        defectPairList =self.GenerateDefectsFromError(7,gateError,tIndex,zIndex,xIndex)
                        probability=self.CNOTErrors[gateError]
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
                    #Second Layer Entangling Gate
                    for gateError in self.CNOTErrors:
                        defectPairList =self.GenerateDefectsFromError(8,gateError,tIndex,zIndex,xIndex)
                        probability=self.CNOTErrors[gateError]
                        if probability>0:
                            AddToNodeGraph(probability,defectPairList)
        for defect1,defect2 in distanceGraph.edges():
            distanceGraph[defect1][defect2]['weight']=-math.log(distanceGraph[defect1][defect2]['weight']/(1-distanceGraph[defect1][defect2]['weight']))
        distanceDict = dict(nx.all_pairs_dijkstra_path_length(distanceGraph))
        return distanceDict
 
    def MatchErrors(self):
        """
        Call after GenerateErrors().
        Generates a weighted graph connecting the error syndromes to each other and the boundaries,
        then stores the minimum weight matching in self.decodedMatches
        """
        evenErrorNodes=[]
        oddErrorNodes=[]
        for point in self.correctMatches.keys():
            if type(point)!=str: # Syndrome node and not a 'node' corresponding to an edge
                if point[1]%2==0 and not (point in evenErrorNodes):
                    evenErrorNodes.append(point)
                elif point[1]%2==1 and not (point in oddErrorNodes):
                    oddErrorNodes.append(point)
        evenErrorGraph=[]
        oddErrorGraph=[]
        evenEdgeNodes=[]
        oddEdgeNodes=[]
        for nodeIndex1 in range(len(evenErrorNodes)):
            for nodeIndex2 in range(nodeIndex1+1,len(evenErrorNodes)):
                weight = self.distanceDict[evenErrorNodes[nodeIndex1]][evenErrorNodes[nodeIndex2]]
                evenErrorGraph.append((evenErrorNodes[nodeIndex1],evenErrorNodes[nodeIndex2],weight))
            minWeight = min([self.distanceDict[evenErrorNodes[nodeIndex1]]["T"],self.distanceDict[evenErrorNodes[nodeIndex1]]["B"]])
            if minWeight==self.distanceDict[evenErrorNodes[nodeIndex1]]["T"]:
                label = ("T",nodeIndex1)
            elif minWeight==self.distanceDict[evenErrorNodes[nodeIndex1]]["B"]:
                label = ("B",nodeIndex1)
            evenErrorGraph.append((evenErrorNodes[nodeIndex1],label,minWeight))
            evenEdgeNodes.append(label)
        for nodeIndex1 in range(len(evenEdgeNodes)):
            for nodeIndex2 in range(nodeIndex1+1,len(evenEdgeNodes)):
                evenErrorGraph.append((evenEdgeNodes[nodeIndex1],evenEdgeNodes[nodeIndex2],0))
        maxWeight = max([abs(graphEdge[2]) for graphEdge in evenErrorGraph]+[0])
        if maxWeight>0:
            scaling = INFTY/(10*maxWeight)
        else:
            scaling = INFTY
        evenErrorGraph = [(x[0],x[1],int(x[2]*scaling)) for x in evenErrorGraph]
        for nodeIndex1 in range(len(oddErrorNodes)):
            for nodeIndex2 in range(nodeIndex1+1,len(oddErrorNodes)):
                weight = self.distanceDict[oddErrorNodes[nodeIndex1]][oddErrorNodes[nodeIndex2]]
                oddErrorGraph.append((oddErrorNodes[nodeIndex1],oddErrorNodes[nodeIndex2],weight))
            minWeight = min([self.distanceDict[oddErrorNodes[nodeIndex1]]["L"],self.distanceDict[oddErrorNodes[nodeIndex1]]["R"]])
            if minWeight==self.distanceDict[oddErrorNodes[nodeIndex1]]["L"]:
                label = ("L",nodeIndex1)
            elif minWeight==self.distanceDict[oddErrorNodes[nodeIndex1]]["R"]:
                label = ("R",nodeIndex1)
            oddErrorGraph.append((oddErrorNodes[nodeIndex1],label,minWeight))
            oddEdgeNodes.append(label)
        for nodeIndex1 in range(len(oddEdgeNodes)):
            for nodeIndex2 in range(nodeIndex1+1,len(oddEdgeNodes)):
                oddErrorGraph.append((oddEdgeNodes[nodeIndex1],oddEdgeNodes[nodeIndex2],0))
        maxWeight = max([abs(graphEdge[2]) for graphEdge in oddErrorGraph]+[0])
        if maxWeight>0:
            scaling = INFTY/(10*maxWeight)
        else:
            scaling = INFTY
        oddErrorGraph = [(x[0],x[1],int(x[2]*scaling)) for x in oddErrorGraph]

        evenMatches = mwpm(evenErrorGraph)
        oddMatches = mwpm(oddErrorGraph)
        self.decodedMatches = evenMatches.union(oddMatches)

    def AttemptCorrectionOfErrors(self):
        """
        Call after MatchErrors.
        Adds the MWPM matches to the actual error graph, to give the error graph after correction
        """
        for defect1,defect2 in self.decodedMatches:
            if type(defect1[0])==str:
                defect1=defect1[0]
            if type(defect2[0])==str:
                defect2=defect2[0]
            if type(defect2)!=str or type(defect1)!=str: 
                self.AddMatchedPair(defect1,defect2)
    
    def FindLogicalErrors(self):
        """
        Call after AttemptCorrectionOfErrors
        Returns the logical errors we made by attempting error correction
        """
        TLogicalError,ZLogicalError,XLogicalError=False,False,False
        if "L" in self.correctMatches.keys():
            if self.correctMatches["L"]=="R":
                ZLogicalError=not ZLogicalError
        if "T" in self.correctMatches.keys():
            if self.correctMatches["T"]=="B":
                XLogicalError=not XLogicalError
        if self.numCrossingsInTime%2==1:
            TLogicalError=True
        return TLogicalError,ZLogicalError,XLogicalError
        
                


"""
dz,dx=10,10
dt=10
eta=10000
p = .005
#S = Syndrome(dz,dx,dt,p,eta,clusterType="RHG")
S = Syndrome(dz,dx,dt,p,eta)
totalT,totalZ,totalX = 0,0,0
for i in range(100):
    S.Clear()
    S.GenerateErrors()
    saveFile=open("saveFile.pk",'wb')
    pickle.dump(S,saveFile)
    saveFile.close()
    #saveFile = open("saveFile.pk",'rb')
    #S=pickle.load(saveFile)
    #saveFile.close()
    S.MatchErrors()
    #fig = plt.figure()
    #ax = fig.add_subplot(1,1,1,projection='3d')
    #S.PlotCluster(ax,plotScaffold=True)
    S.AttemptCorrectionOfErrors()
    T,Z,X = S.FindLogicalErrors()
    if Z:
        totalZ+=1
    if X:
        totalX+=1
    if T:
        totalT+=1
    plt.show()
print(totalT,totalZ,totalX)
"""
