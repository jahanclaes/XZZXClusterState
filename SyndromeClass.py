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

    def __init__(self,p,eta):
        """
        Initializes a Syndrome object
        dx and dz are the dimensions of the code
        numMeasurements is the number of rounds of stabilizer measurements we are simulating
        """
        self.dx=dx #X-distance of the cluster
        self.dz=dz #Z-distance of the cluster
        self.dt = dt #Time-distance of the cluster
        self.p = p #Probability of Z errors
        self.eta = eta #Probability of X errors is p/eta
        self.correctMatches = {}
        self.decodedMatches = {}

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
                    #Ancilla qubit errors
                    error = random.choices(["I","X","Y","Z"],weights = [1-self.p*(1-2/self.eta),self.p/self.eta,self.p/self.eta,self.p])
                    if error == "Y" or error=="Z": #Generate a defect pair
                        if tIndex==0:
                            self.AddMatchedPair("B",(tIndex,zIndex,xIndex)) #"B" is the node corresponding to the bottom edge of the cluster
                        elif tIndex==self.dt-1:
                            self.AddMatchedPair("T",(tIndex,zIndex,xIndex)) #"T" is the node corresponding to the top edge of the cluster
                        else:
                            self.AddMatchedPair((tIndex-1,zIndex,xIndex),(tIndex,zIndex,xIndex))
                for xIndex in range(self.dx-(zIndex)%2):
                    #Data qubit errors
                    error = random.choices(["I","X","Y","Z"],weights = [1-self.p*(1-2/self.eta),self.p/self.eta,self.p/self.eta,self.p]) #Error on lower qubit pair
                    
                    error = random.choices(["I","X","Y","Z"],weights = [1-self.p*(1-2/self.eta),self.p/self.eta,self.p/self.eta,self.p]) #Error on upper qubit pair
             
                    
                        
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
                


        
"""
dz,dx=2,3
numMeasurements=5
eta=1
p = .005
S = Syndrome(dz,dx,numMeasurements,p,eta)
S.GenerateDistanceDict()
for i in range(100):
    S = Syndrome(dz,dx,numMeasurements,p,eta)
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
""" 
