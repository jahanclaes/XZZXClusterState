import numpy as np
import matplotlib
matplotlib.use('TkAgg')
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

    def __init__(self,dt,dz,dx,p,eta,XZZX=True):
        """
        Initializes a Syndrome object
        dx and dz are the dimensions of the code
        """
        self.dx=dx #X-distance of the cluster
        self.dz=dz #Z-distance of the cluster
        self.dt = dt #Time-distance of the cluster
        self.p = p #Probability of Z errors
        self.eta = eta #Probability of X errors is p/eta
        self.XZZX = XZZX
        self.correctMatches = {}
        self.errorChains = []
        self.decodedMatches = []

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
                point2 = (point1[0],point1[1],self.dx-1)
            elif point2[0]=="B" and type(point1)!=str:
                point2 = (point1[0],point1[1],-1)
            elif point2[0]=="L" and type(point1)!=str:
                point2 = (point1[0],-1,point1[2])
            elif point2[0]=="R" and type(point1)!=str:
                point2 = (point1[0],2*self.dz-1,point1[2])
            if type(point1)!=str:
                if point1[1]%2==0:
                    plt.plot([point1[1]/2,point2[1]/2],[point1[2]+((point1[1]+1)%2)/2,point2[2]+((point2[1]+1)%2)/2],[point1[0]+.5,point2[0]+.5],color=matchColor,lw=errorWeight,ls='dotted')
                else:
                    plt.plot([point1[1]/2,point2[1]/2],[point1[2]+((point1[1]+1)%2)/2,point2[2]+((point2[1]+1)%2)/2],[point1[0]+1,point2[0]+1],color=matchColor,lw=errorWeight,ls='dotted')
                
                
                    

    def AddMatchedPair(self,point1,point2):
        if point1!=point2:
            self.errorChains.append((point1,point2))
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

    def GenerateErrors(self):
        """
        Generate defects according to some error model.
        Right now, we are assuming phenomenological biased noise occurs after initializing the cluster state.
        If we instead used circuit level noise, would lead to correlations but shouldn't increase the bias
        Populates self.correctMatches
        """
        for tIndex in range(self.dt):
            for zIndex in range(2*self.dz-1):
                for xIndex in range(self.dx-1+zIndex%2):
                    #Error on ancilla qubit
                    error = random.choices(["I","X","Y","Z"],weights = [1-self.p*(1-2/self.eta),self.p/self.eta,self.p/self.eta,self.p])[0]
                    if error == "Z" or error=="Y": #Generate a defect pair
                        self.AddMatchedPair((tIndex-1,zIndex,xIndex),(tIndex,zIndex,xIndex))
                for xIndex in range(self.dx-(zIndex)%2):
                    #Error on lower data qubit
                    error = random.choices(["I","X","Y","Z"],weights = [1-self.p*(1-2/self.eta),self.p/self.eta,self.p/self.eta,self.p])[0]
                    if zIndex%2==0 and (error=="Z" or error=="Y"):
                        if zIndex==0:
                            self.AddMatchedPair("L",(tIndex-1,zIndex+1,xIndex)) #L is the node corresponding to the left of the cluster
                        elif zIndex==2*self.dz-2:
                            self.AddMatchedPair("R",(tIndex-1,zIndex-1,xIndex)) #R is the node corresponding to the right of the cluster
                        else:
                            self.AddMatchedPair((tIndex-1,zIndex-1,xIndex),(tIndex-1,zIndex+1,xIndex))
                    elif zIndex%2==1 and (((error=="X" or error=="Y") and self.XZZX) or ((error=="Z" or error=="Y") and not self.XZZX)):
                        self.AddMatchedPair((tIndex-1,zIndex,xIndex),(tIndex-1,zIndex,xIndex+1))
                    #Error on upper data qubit
                    error = random.choices(["I","X","Y","Z"],weights = [1-self.p*(1-2/self.eta),self.p/self.eta,self.p/self.eta,self.p])[0]
                    if zIndex%2==0 and (((error=="X" or error=="Y") and self.XZZX) or ((error=="Z" or error=="Y") and not self.XZZX)):
                        if xIndex==0:
                            self.AddMatchedPair("B",(tIndex,zIndex,xIndex)) #B is the node corresponding to the bottom of the cluster
                        elif xIndex==self.dx-(zIndex)%2-1:
                            self.AddMatchedPair("T",(tIndex,zIndex,xIndex-1)) #T is the node corresponding to the top of the cluster
                        else:
                            self.AddMatchedPair((tIndex,zIndex,xIndex),(tIndex,zIndex,xIndex-1))
                    elif zIndex%2==1 and (error=="Z" or error=="Y"):
                        self.AddMatchedPair((tIndex,zIndex-1,xIndex),(tIndex,zIndex+1,xIndex))
                                           
    def ComputeWeight(self,node1,node2):
        p=self.p
        eta=self.eta
        if node2=="T" and self.XZZX:
            weight = -math.log(2*p/eta/(1-2*p/eta))*(self.dx-1-node1[2])
        elif node2=="T" and not self.XZZX:
            weight = -math.log((p+p/eta)/(1-p-p/eta))*(self.dx-1-node1[2])
        elif node2=="B" and self.XZZX:
            weight = -math.log(2*p/eta/(1-2*p/eta))*(node1[2]+1)
        elif node2=="B" and not self.XZZX:
            weight = -math.log((p+p/eta)/(1-p-p/eta))*(node1[2]+1)
        elif node2=="L":
            weight = -math.log((p+p/eta)/(1-p-p/eta))*(node1[1]+1)/2
        elif node2=="R":
            weight = -math.log((p+p/eta)/(1-p-p/eta))*(2*self.dz-node1[1]-1)/2
        else:
            zDistance =abs(node2[1]-node1[1])/2
            tDistance =abs(node2[0]-node1[0])
            xDistance =abs(node2[2]-node1[2]) 
            if self.XZZX:
                weight = -(math.log((p+p/eta)/(1-p-p/eta))*(zDistance+tDistance)+math.log(2*p/eta/(1-2*p/eta))*xDistance)
            else:
                weight = -(math.log((p+p/eta)/(1-p-p/eta))*(zDistance+tDistance+xDistance))
        return weight

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
                weight = self.ComputeWeight(evenErrorNodes[nodeIndex1],evenErrorNodes[nodeIndex2])
                evenErrorGraph.append((evenErrorNodes[nodeIndex1],evenErrorNodes[nodeIndex2],weight))
            edgeWeights = [self.ComputeWeight(evenErrorNodes[nodeIndex1],edge) for edge in ["T","B"]]
            minWeight = min(edgeWeights)
            label = ["T","B"][edgeWeights.index(minWeight)]+str(nodeIndex1)
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
                weight = self.ComputeWeight(oddErrorNodes[nodeIndex1],oddErrorNodes[nodeIndex2])
                oddErrorGraph.append((oddErrorNodes[nodeIndex1],oddErrorNodes[nodeIndex2],weight))
            edgeWeights = [self.ComputeWeight(oddErrorNodes[nodeIndex1],edge) for edge in ["L","R"]]
            minWeight = min(edgeWeights)
            label = ["L","R"][edgeWeights.index(minWeight)]+str(nodeIndex1)
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
        ZLogicalError,XLogicallError=False,False
        if "L" in self.correctMatches.keys():
            if self.correctMatches["L"]=="R":
                ZLogicalError=not ZLogicalError
        if "T" in self.correctMatches.keys():
            if self.correctMatches["T"]=="B":
                XLogicallError=not XLogicallError
        return ZLogicalError,XLogicallError
                


dz,dx,dt=3,3,3
eta=100
p = .05
S = Syndrome(dt,dz,dx,p,eta)
fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')
S.PlotCluster(ax)
plt.show()
""" 
for i in range(10):
    S = Syndrome(dt,dz,dx,p,eta)
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
    S.PlotCluster(ax,plotScaffold=False)
    S.AttemptCorrectionOfErrors()
    print("Z,X",S.FindLogicalErrors())
    plt.show()
"""
