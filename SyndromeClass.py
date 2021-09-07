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
Overall measurement pattern on a single stabilizer is X-Z-Z-X:
  Z      2
 / \    / \ 
X   X  1   4 
 \ /    \ /
  Z      3  
"""

class Syndrome:

    def __init__(self,dz,dx,dt,p,eta):
        """
        Initializes a Syndrome object
        dx and dz are the dimensions of the code
        dt is the number of rounds of stabilizer measurements we are simulating
        """
        self.dx=dx
        self.dz=dz
        self.dt = dt
        self.freePlaquettesStart = [(2*dz-1,xIndex) for xIndex in range(dx)]
        self.freePlaquettesEnd = [(2*dz-1,xIndex) for xIndex in range(dx)]
        self.correctMatches = {}
        self.decodedMatches = {}
        pIdle=p
        pPrep=p
        pMeasure=p
        pCNOT=p
        pCZ=p
        singleQubitErrors = ["I","Z","X","Y"]
        twoQubitErrors =["II","IZ","ZI","ZZ","IX","IY","XI","XX","XY","XZ","YI","YX","YY","YZ","ZX","ZY"] 
        weightsSPAM = [1-pMeasure-pPrep,pMeasure+pPrep,0,0]
        weightsIdle = [0]*4 
        weightsCNOT = [0]*16
        weightsCZ = [0]*16
        self.SPAMErrors = {singleQubitErrors[i]:weightsSPAM[i] for i in range(4)}
        self.idleErrors = {singleQubitErrors[i]:weightsIdle[i] for i in range(4)}
        self.CNOTErrors = {twoQubitErrors[i]:weightsCNOT[i] for i in range(16)}
        self.CZErrors = {twoQubitErrors[i]:weightsCZ[i] for i in range(16)}
        try:
            saveFile = open("distanceDict_"+str(self.dt)+"_"+str(self.dz)+"_"+str(self.dx)+"_"+str(self.CNOTErrors["ZZ"])+".pk",'rb')
            self.distanceDict = pickle.load(saveFile)
            saveFile.close()
            print("Loaded")
        except:
            self.distanceDict = self.GenerateDistanceDict()
 

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
        for point1,point2 in self.errorChains:
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
            if type(point1)==str:
                point1,point2=point2,point1
            if point2[0]=="T" and type(point1)!=str:

    def AddMatchedPair(self,point1,point2):
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

    def GenerateDefectsFromError(self,errorIndex,errorString,tIndex,zIndex,xIndex):
        defectPairList=[]
        #Ancilla Preparation/Measurement errors
        if errorIndex==0:
            if errorString == "Z" or errorString=="Y":
                defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex,xIndex)))
        #Idle Data Qubit Errors
        elif errorIndex==1:
            if zIndex==1:
                if errorString == "Z" or errorString=="Y":
                    defectPairList.append(((tIndex+1,zIndex,xIndex),(tIndex+1,zIndex-2,xIndex)))
                if errorString == "X" or errorString=="Y":
                    defectPairList.append(((tIndex+1,zIndex-1,xIndex-zIndex%2),(tIndex+1,zIndex-1,xIndex-zIndex%2+1)))
            if zIndex==4*self.dz-3:
                if errorString == "Z" or errorString=="Y":
                    defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex+2,xIndex)))
                if errorString == "X" or errorString=="Y":
                    defectPairList.append(((tIndex,zIndex+1,xIndex-zIndex%2),(tIndex,zIndex+1,xIndex-zIndex%2+1)))
            if xIndex==0 and zIndex%2==0:
                if errorString == "Z" or errorString=="Y":
                    defectPairList.append(((tIndex,zIndex-1,xIndex-zIndex%2),(tIndex+1,zIndex+1,xIndex-zIndex%2)))
                if errorString == "X" or errorString=="Y":
                    defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex,xIndex-1)))
            if xIndex==self.dx-2 and zIndex%2==0:
                if errorString == "Z" or errorString=="Y":
                    defectPairList.append(((tIndex,zIndex-1,xIndex-zIndex%2+1),(tIndex+1,zIndex+1,xIndex-zIndex%2+1)))
                if errorString == "X" or errorString=="Y":
                    defectPairList.append(((tIndex+1,zIndex,xIndex),(tIndex+1,zIndex,xIndex+1)))
        #Idle Ancilla Qubit Errors
        elif errorIndex==2:
            if zIndex==0:
                pass # Can just prepare the state one timestep later
            if zIndex==4*self.dz-2:
                pass # Can just measure the state one timestep earlier
            if xIndex==0 and zIndex%2==1:
                if errorString == "X" or errorString=="Y":
                    defectPairList.append(((tIndex+1,zIndex+1,xIndex-zIndex%2+1),(tIndex+1,zIndex+1,xIndex-zIndex%2)))
                if errorString == "Z" or errorString=="Y":
                    defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex,xIndex)))
            if xIndex==self.dx-1 and zIndex%2==1:
                if errorString == "X" or errorString=="Y":
                    defectPairList.append(((tIndex,zIndex-1,xIndex-zIndex%2),(tIndex+1,zIndex+1,xIndex-zIndex%2+1)))
                if errorString == "Z" or errorString=="Y":
                    defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex,xIndex)))
        #First 2-qubit (CNOT) Error
        elif errorIndex==3:
            if zIndex>0:
                if errorString[0]=="X" or errorString[0]=="Y": #Data (Target) qubit error
                    defectPairList.append(((tIndex,zIndex-1,xIndex+1-zIndex%2),(tIndex,zIndex-1,xIndex-zIndex%2)))
                if errorString[0]=="Z" or errorString[0]=="Y": #Data (Target) qubit error
                    defectPairList.append(((tIndex,zIndex-2,xIndex),(tIndex+1,zIndex,xIndex)))
                if errorString[1]=="X" or errorString[1]=="Y": #Ancilla (Control) qubit error
                    defectPairList.append(((tIndex,zIndex-1,xIndex),(tIndex,zIndex-1,xIndex+1-2*(zIndex%2))))
                if errorString[1]=="Z" or errorString[1]=="Y": #Ancilla (Control) qubit error
                    defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex,xIndex)))
        #Second 2-qubit (CZ) Error
        elif errorIndex==4:
            if xIndex<self.dx-1:
                if errorString[0]=="X" or errorString[0]=="Y": #Data (Target) qubit error
                    defectPairList.append(((tIndex,zIndex,xIndex+1),(tIndex+1,zIndex,xIndex)))
                if errorString[0]=="Z" or errorString[0]=="Y": #Data (Target) qubit error
                    defectPairList.append(((tIndex,zIndex-1,xIndex-zIndex%2+1),(tIndex+1,zIndex+1,xIndex-zIndex%2+1)))
                if errorString[1]=="X" or errorString[1]=="Y": #Ancilla (Control) qubit error
                    defectPairList.append(((tIndex,zIndex-1,xIndex-zIndex%2),(tIndex+1,zIndex+1,xIndex-zIndex%2+1)))
                if errorString[1]=="Z" or errorString[1]=="Y": #Ancilla (Control) qubit error
                    defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex,xIndex)))
        #Third 2-qubit (CZ) Error
        elif errorIndex==5:
            if xIndex>zIndex%2-1:
                if errorString[0]=="X" or errorString[0]=="Y": #Data (Target) qubit error
                    defectPairList.append(((tIndex+1,zIndex,xIndex),(tIndex+1,zIndex,xIndex-1)))
                if errorString[0]=="Z" or errorString[0]=="Y": #Data (Target) qubit error
                    defectPairList.append(((tIndex,zIndex-1,xIndex-zIndex%2),(tIndex+1,zIndex+1,xIndex-zIndex%2)))
                if errorString[1]=="X" or errorString[1]=="Y": #Ancilla (Control) qubit error
                    defectPairList.append(((tIndex+1,zIndex+1,xIndex-zIndex%2+1),(tIndex+1,zIndex+1,xIndex-zIndex%2)))
                if errorString[1]=="Z" or errorString[1]=="Y": #Ancilla (Control) qubit error
                    defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex,xIndex)))
        #Fourth 2-qubit (CNOT) Error
        elif errorIndex==6:
            if zIndex<4*self.dz-2:
                if errorString[0]=="X" or errorString[0]=="Y": #Data (Target) qubit error
                    defectPairList.append(((tIndex+1,zIndex+1,xIndex-zIndex%2),(tIndex+1,zIndex+1,xIndex-zIndex%2+1)))
                if errorString[0]=="Z" or errorString[0]=="Y": #Data (Target) qubit error
                    defectPairList.append(((tIndex+1,zIndex,xIndex),(tIndex+1,zIndex+2,xIndex)))
                if errorString[1]=="X" or errorString[1]=="Y": #Ancilla (Control) qubit error
                    pass
                if errorString[1]=="Z" or errorString[1]=="Y": #Ancilla (Control) qubit error
                    defectPairList.append(((tIndex,zIndex,xIndex),(tIndex+1,zIndex,xIndex)))
        for index in range(len(defectPairList)):
            point1,point2 = defectPairList[index]
            if point1[1]==-1:
                point1="L"
            elif point1[1]==4*self.dz-1:
                point1="R"
            elif point1[2]==-1:
                point1="B"
            elif point1[2]==self.dx-1 and point1[1]%2==0:
                point1="T"
            elif point1[0]==0 and (point1[1],point1[2]) in self.freePlaquettesStart:
                point1="S"
            elif point1[0]==self.dt-1 and (point1[1],point1[2]) in self.freePlaquettesEnd:
                point1="E"
            if point2[1]==-1:
                point2="L"
            elif point2[1]==4*self.dz-1:
                point2="R"
            elif point2[2]==-1:
                point2="B"
            elif point2[2]==self.dx-1 and point2[1]%2==0:
                point2="T"
            elif point2[0]==0 and (point2[1],point2[2]) in self.freePlaquettesStart:
                point2="S"
            elif point2[0]==self.dt-1 and (point2[1],point2[2]) in self.freePlaquettesEnd:
                point2="E"
            defectPairList[index]=(point1,point2)
        return defectPairList

    def GenerateErrors(self):
        """
        Generate defects according to some error model.
        Populates self.correctMatches
        """
        for tIndex in range(self.dt-1):
            for zIndex in range(4*self.dz-1):
                for xIndex in range(self.dx-1+zIndex%2):
                    #Ancilla Preparation/Measurement errors
                    defectPairList = []
                    ancillaError = random.choices(list(self.SPAMErrors.keys()),weights=[self.SPAMErrors[key] for key in self.SPAMErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(0,ancillaError,tIndex,zIndex,xIndex)
                    #Idle Data Qubit Errors
                    idleError = random.choices(list(self.idleErrors.keys()),weights = [self.idleErrors[key] for key in self.idleErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(1,idleError,tIndex,zIndex,xIndex)
                    #Idle Ancilla Qubit Errors
                    idleError = random.choices(list(self.idleErrors.keys()),weights = [self.idleErrors[key] for key in self.idleErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(2,idleError,tIndex,zIndex,xIndex)
                    #First 2-qubit (CNOT) Error
                    gateError = random.choices(list(self.CNOTErrors.keys()),weights = [self.CNOTErrors[key] for key in self.CNOTErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(3,gateError,tIndex,zIndex,xIndex)
                    #Second 2-qubit (CZ) Error
                    gateError = random.choices(list(self.CZErrors.keys()),weights = [self.CZErrors[key] for key in self.CZErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(4,gateError,tIndex,zIndex,xIndex)
                    #Third 2-qubit (CZ) Error
                    gateError = random.choices(list(self.CZErrors.keys()),weights = [self.CZErrors[key] for key in self.CZErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(5,gateError,tIndex,zIndex,xIndex)
                    #Fourth 2-qubit (CNOT) Error
                    gateError = random.choices(list(self.CNOTErrors.keys()),weights = [self.CNOTErrors[key] for key in self.CNOTErrors.keys()])[0]
                    defectPairList = defectPairList+self.GenerateDefectsFromError(6,gateError,tIndex,zIndex,xIndex)
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

        for tIndex in range(self.dt-1):
            for zIndex in range(4*self.dz-1):
                for xIndex in range(self.dx-1+zIndex%2):
                    #Ancilla Preparation/Measurement errors
                    for ancillaError in self.SPAMErrors:
                        defectPairList = self.GenerateDefectsFromError(0,ancillaError,tIndex,zIndex,xIndex)
                        if self.SPAMErrors[ancillaError]>0:
                            probability = self.SPAMErrors[ancillaError]
                            AddToNodeGraph(probability,defectPairList)
                    #Idle Data Qubit Errors
                    for idleError in self.idleErrors:
                        defectPairList = self.GenerateDefectsFromError(1,idleError,tIndex,zIndex,xIndex)
                        probability = (self.idleErrors[idleError])
                        AddToNodeGraph(probability,defectPairList)
                    #Idle Ancilla Qubit Errors
                    for idleError in self.idleErrors:
                        defectPairList = self.GenerateDefectsFromError(2,idleError,tIndex,zIndex,xIndex)
                        probability = (self.idleErrors[idleError])
                        AddToNodeGraph(probability,defectPairList)
                    #First 2-qubit (CNOT) Error
                    for gateError in self.CNOTErrors:
                        defectPairList = self.GenerateDefectsFromError(3,gateError,tIndex,zIndex,xIndex)
                        probability = (self.CNOTErrors[gateError])
                        AddToNodeGraph(probability,defectPairList)
                    #Second 2-qubit (CZ) Error
                    for gateError in self.CZErrors:
                        defectPairList = self.GenerateDefectsFromError(4,gateError,tIndex,zIndex,xIndex)
                        probability = (self.CZErrors[gateError])
                        AddToNodeGraph(probability,defectPairList)
                    #Third 2-qubit (CZ) Error
                    for gateError in self.CZErrors:
                        defectPairList = self.GenerateDefectsFromError(5,gateError,tIndex,zIndex,xIndex)
                        probability = (self.CZErrors[gateError])
                        AddToNodeGraph(probability,defectPairList)
                    #Fourth 2-qubit (CNOT) Error
                    for gateError in self.CNOTErrors:
                        defectPairList = self.GenerateDefectsFromError(6,gateError,tIndex,zIndex,xIndex)
                        probability = self.CNOTErrors[gateError]
                        AddToNodeGraph(probability,defectPairList)
        for defect1,defect2 in distanceGraph.edges():
            distanceGraph[defect1][defect2]['weight']=-math.log(distanceGraph[defect1][defect2]['weight'])
        distanceDict = dict(nx.all_pairs_dijkstra_path_length(distanceGraph))
        saveFile = open("distanceDict_"+str(self.dt)+"_"+str(self.dz)+"_"+str(self.dx)+"_"+str(self.CNOTErrors["ZZ"])+".pk",'wb')
        pickle.dump(distanceDict,saveFile)
        saveFile.close()
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
            if "S" in self.distanceDict[evenErrorNodes[nodeIndex1]]:
                minWeight = min([self.distanceDict[evenErrorNodes[nodeIndex1]]["T"],self.distanceDict[evenErrorNodes[nodeIndex1]]["B"],self.distanceDict[evenErrorNodes[nodeIndex1]]["S"],self.distanceDict[evenErrorNodes[nodeIndex1]]["E"]])
                if minWeight==self.distanceDict[evenErrorNodes[nodeIndex1]]["T"]:
                    label = ("T",nodeIndex1)
                elif minWeight==self.distanceDict[evenErrorNodes[nodeIndex1]]["B"]:
                    label = ("B",nodeIndex1)
                elif minWeight==self.distanceDict[evenErrorNodes[nodeIndex1]]["S"]:
                    label = ("S",nodeIndex1)
                elif minWeight==self.distanceDict[evenErrorNodes[nodeIndex1]]["E"]:
                    label = ("E",nodeIndex1)
            else:
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
            minWeight = min([self.distanceDict[oddErrorNodes[nodeIndex1]]["L"],self.distanceDict[oddErrorNodes[nodeIndex1]]["R"],self.distanceDict[oddErrorNodes[nodeIndex1]]["S"],self.distanceDict[oddErrorNodes[nodeIndex1]]["E"]])
            if minWeight==self.distanceDict[oddErrorNodes[nodeIndex1]]["L"]:
                label = ("L",nodeIndex1)
            elif minWeight==self.distanceDict[oddErrorNodes[nodeIndex1]]["R"]:
                label = ("R",nodeIndex1)
            elif minWeight==self.distanceDict[oddErrorNodes[nodeIndex1]]["S"]:
                label = ("S",nodeIndex1)
            elif minWeight==self.distanceDict[oddErrorNodes[nodeIndex1]]["E"]:
                label = ("E",nodeIndex1)
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
        ZLogicalError,Z1LogicalError,Z2LogicalError,XLogicallError,TimeLogicalError=False,False,False,False,False
        if "S" in self.correctMatches.keys():
            if self.correctMatches["S"]=="L":
                Z1LogicalError=not Z1LogicalError
            if self.correctMatches["S"]=="R":
                Z2LogicalError=not Z2LogicalError
            if self.correctMatches["S"]=="E":
                TimeLogicalError=not TimeLogicalError
        if "E" in self.correctMatches.keys():
            if self.correctMatches["E"]=="L":
                Z1LogicalError=not Z1LogicalError
            if self.correctMatches["E"]=="R":
                Z2LogicalError=not Z2LogicalError
        if "L" in self.correctMatches.keys():
            if self.correctMatches["L"]=="R":
                ZLogicalError=not ZLogicalError
        if "T" in self.correctMatches.keys():
            if self.correctMatches["T"]=="B":
                XLogicallError=not XLogicallError
        return ZLogicalError,Z1LogicalError,Z2LogicalError,XLogicallError,TimeLogicalError
                


        
dz,dx=2,3
dt=5
eta=1
p = .005
S = Syndrome(dz,dx,dt,p,eta)
S.GenerateDistanceDict()
for i in range(100):
    S = Syndrome(dz,dx,dt,p,eta)
    S.GenerateErrors()
    saveFile=open("saveFile.pk",'wb')
    pickle.dump(S,saveFile)
    saveFile.close()
    #saveFile = open("saveFile.pk",'rb')
    #S=pickle.load(saveFile)
    #saveFile.close()
    S.MatchErrors()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection='3d')
    S.Plot3D(ax,parity=0)
    S.Plot3D(ax,parity=1)
    S.AttemptCorrectionOfErrors()
    print(S.correctMatches)
    print("Z,Z1,Z2,X,T",S.FindLogicalErrors())
    plt.show()
