###########################################
#             Thai T. Pham                #
#         Stanford University             #
###########################################

# * Warning: The codes are not optimized for speed.

# Exact and Approximate betweenness centrality

import snap
import Queue
from copy import deepcopy
import random

def ComputeBetweennessCentrality(G, sNode, B, type, k):
	s = sNode.GetId() 
	sigma = {} # number of shortest paths between sNode and other nodes
	distance = {}     # distance between each node and sNode
	ParentSet = {}
	ChildrenSet = {}

	for node in G.Nodes():
		i = node.GetId()
		if i == s:
			sigma[i] = 1
			distance[i] = 0
		else:
			sigma[i] = 0
			distance[i] = -1
		ParentSet[i] = []
		ChildrenSet[i] = []

	Q = Queue.Queue()
	S = [] # set of nodes through BFS from sNode
	Q.put(s)

	# compute ParentSet, ChildrenSet, and sigma recursively through BFS
	while not Q.empty():
		u = Q.get()
		S.append(u)
		for v in G.GetNI(u).GetOutEdges():
			if distance[v] < 0:
				Q.put(v)
				distance[v] = distance[u] + 1
			if distance[v] == distance[u] + 1:
				sigma[v] = sigma[v] + sigma[u]
				ParentSet[v].append(u)
				ChildrenSet[u].append(v)

	# delta: dependency scores of sNode on each edge
	delta = deepcopy(B) # initialize delta

	# compute delta and B recursively
	while len(S) > 0:
		w = S.pop()
		for v in ParentSet[w]:
			vw = (min(v, w), max(v, w))
			if len(ChildrenSet[w]) == 0: # w is a leaf
				delta[vw] = sigma[v] * 1.0 / sigma[w]		
			else:
				delta[vw] = 1
				for x in ChildrenSet[w]:
					wx = (min(w, x), max(w, x))
					delta[vw] = delta[vw] + delta[wx]
				delta[vw] = sigma[v] * 1.0 / sigma[w] * delta[vw]

			if type == 1:
				B[vw] = B[vw] + delta[vw]
			else:
				if B[vw] <= c * G.GetNodes():
					B[vw] = B[vw] + delta[vw]
					k[vw] = k[vw] + 1



# Simulation
n = 1000
c = 5

Rnd = snap.TRnd()
G = snap.GenPrefAttach(n, 4, Rnd)
B1 = {}  # Betweenness Centrality for algorithm 1
B2 = {}	 # Betweenness Centrality for algorithm 2
k = {}   # Number of samples for each edge
edges = []
for edge in G.Edges():
	edges.append((edge.GetSrcNId(), edge.GetDstNId()))
for i in xrange(len(edges)):
	B1[edges[i]] = 0 
	B2[edges[i]] = 0 
	k[edges[i]] = 0  

# compute betweenness centrality for algorithm 1 - exact betweenness centrality
for sNode in G.Nodes():
	ComputeBetweennessCentrality(G, sNode, B1, 1, k)

# compute betweennness centrality for algorithm 2 - approximate betweenness centrality
numSamples = n / 10
nodeIds = [node.GetId() for node in G.Nodes()]
for i in xrange(numSamples):
	sourceNId = random.choice(nodeIds)
	s = G.GetNI(sourceNId)
	nodeIds.remove(sourceNId)
	ComputeBetweennessCentrality(G, s, B2, 2, k)

# report result
BC1 = []
BC2 = []
for i in xrange(len(edges)):
	BC1.append(B1[edges[i]])
	BC2.append(n * 1.0 / k[edges[i]] * B2[edges[i]])

BC1.sort(reverse = True)
BC2.sort(reverse = True)
xVals = [x for x in xrange(1, len(BC1) + 1)]


