###########################################
#             Thai T. Pham                #
#         Stanford University             #
###########################################

# * Warning: The code is not optimized for speed.

import snap
import random
import numpy
import math

##################################################
########### Signed Networks Over Time ############
##################################################

# (I) Create a complete network on 10 nodes.
# (II) For each edge, choose a sign (+,−) at random (p = 1/2).
# (III) Run the dynamic process described below for 1,000,000 iterations:
	# (i) Pick a triad at random.
	# (ii) If it’s balanced, do nothing.
	# (iii) Otherwise, choose one of the edges uniformly at random and flip its sign

n = 10
numIter = 1000000
numTime = 100

def generateGplus():
	G = snap.GenFull(snap.PUNGraph, n)
	Gsigns = [[1 for x in xrange(n)] for x in xrange(n)]

	for i in xrange(n - 1):
		for j in xrange(i + 1, n):
			randomValue = random.random()
			if randomValue >= 0.5:
				Gsigns[i][j] = -1
				Gsigns[j][i] = -1

	return (G, Gsigns)

def isBalanced(G, Gsigns):
	nodeSet1 = []
	nodeSet2 = []

	nodeSet1.append(0)
	for i in xrange(1, n):
		if Gsigns[0][i] == 1:
			nodeSet1.append(i)
		else:
			nodeSet2.append(i)

	for x in xrange(len(nodeSet1) - 1):
		for y in xrange(i + 1, len(nodeSet1)):
			i = nodeSet1[x]
			j = nodeSet1[y]
			if Gsigns[i][j] == -1:
				return False

	for x in xrange(len(nodeSet2) - 1):
		for y in xrange(i + 1, len(nodeSet2)):
			i = nodeSet2[x]
			j = nodeSet2[y]
			if Gsigns[i][j] == -1:
				return False

	for x in xrange(len(nodeSet1)):
		for y in xrange(len(nodeSet2)):
			i = nodeSet1[x]
			j = nodeSet2[y]
			if Gsigns[i][j] == 1:
				return False

	return True

def dynamicProcess(G, Gsigns):
	count = 0
	if count < 1:
		i = random.randint(0, n - 1)
		j = random.randint(0, n - 1)
		k = random.randint(0, n - 1)

		if i != j and j != k and k != i:
			count += 1

	if (Gsigns[i][j] + Gsigns[j][k] + Gsigns[k][i] == -3) or \
		(Gsigns[i][j] + Gsigns[j][k] + Gsigns[k][i] == 1):
		randomNum = random.random()
		if randomNum < 1.0 / 3:
			Gsigns[i][j] = - Gsigns[i][j]
			Gsigns[j][i] = Gsigns[i][j]
		elif randomNum < 2.0 / 3:
			Gsigns[j][k] = - Gsigns[j][k]
			Gsigns[k][j] = Gsigns[j][k]
		else:
			Gsigns[k][i] = -Gsigns[k][i]
			Gsigns[i][k] = Gsigns[k][i]

	return (G, Gsigns)

#########################################
count = 0
for i in xrange(numTime):
	(G, Gsigns) = generateGplus()
	for j in xrange(numIter):
		if not isBalanced(G, Gsigns):
			(G, Gsigns) = dynamicProcess(G, Gsigns)

	if isBalanced(G, Gsigns):
		count += 1

print count * 1.0 / numTime # result: 1.0

