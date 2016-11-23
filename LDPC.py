###############################
#         Thai T. Pham        #
#      Stanford University    #
###############################

# NOTICE: The code is not optimized for speed #

from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

# This is an algorithms for reliable communication in the presence of noise. 
# We focus on error correcting codes based on highly sparse, low density parity 
# check (LDPC) matrices, and use the sumproduct variant of the loopy belief 
# propagation (BP) algorithm to estimate partially corrupted message
# bits. For background information on LDPC codes, see Chap. 47 of 
# MacKays Information Theory, Inference, and Learning Algorithms, 
# which is freely available online: http://www.inference.phy.cam.ac.uk/mackay/itila/.

# We consider rate 1/2 error correcting codes, which encode N message bits using a 
# 2N-bit codeword. LDPC codes are specified by a N X 2N binary parity check matrix 
# H, whose columns correspond to codeword bits, and rows to parity check constraints. 
# We define Hij = 1 if parity check i depends on codeword bit j, and Hij = 0 otherwise. 
# Valid codewords are those for which the sum of the bits connected to each parity
# check, as indicated by H, equals zero in modulo-2 addition (i.e., the number of 
# active bits must be even). Equivalently, the modulo-2 product of the parity check 
# matrix with the 2N-bit codeword vector must equal a N-bit vector of zeros. As 
# illustrated in Fig. 1, we can visualize these parity check constraints via
# a corresponding factor graph. The parity check matrix H can then be thought of as 
# an adjacency matrix, where rows correspond to factor (parity) nodes, columns to 
# variable (codeword bit) nodes, and ones to edges linking factors to variables.

# For the LDPC codes we consider, we also define a corresponding 2N X N generator 
# matrix G. To encode an N-bit message vector we would like to transmit, we take 
# the modulo-2 matrix product of the generator matrix with the message. The generator 
# matrix has been constructed (via linear algebra over the finite field GF(2)) such 
# that this product always produces a valid 2N-bit codeword. Geometrically, its 
# columns are chosen to span the null space of H. We use a systematic encoding, in 
# which the first N codeword bits are simply copies of the message bits.

maxIter = 50

class ldpc_graph:
    def __init__(self, H, noisy_output, eps):
        self.H = H
        self.noisy_output = noisy_output
        self.eps = eps

        self.factors = {}
        self.nodes = {}

        self.unary_factors = {}
        self.mar_props = {}
        self.beliefs = []

    # initialize graph data
    def create_ldpc_graph(self):
        H = self.H
        noisy_output = self.noisy_output

        N = H.shape[0]

        for j in xrange(len(noisy_output)):
            if noisy_output[j] == 0:
                self.unary_factors[j] = np.array([1-self.eps, self.eps])
            else:
                self.unary_factors[j] = np.array([self.eps, 1-self.eps])

        for i in xrange(N):
            self.factors[i] = {}
            val = np.dot(H[i, :], noisy_output) % 2
            for j in xrange(2*N):
                if H[i, j] == 1:
                    if val == 0:
                        self.factors[i][j] = np.array([1-self.eps, self.eps])
                    else:
                        self.factors[i][j] = np.array([self.eps, 1-self.eps])

        for j in xrange(2*N):
            self.nodes[j] = {}
            for i in xrange(N):
                if H[i, j] == 1:
                    self.nodes[j][i] = self.unary_factors[j]

    def get_beliefs(self):
        return self.beliefs

    # main function - loopy belief propagation
    def loopy_belief_prop(self, maxIter):
        nodes = self.nodes
        factors = self.factors
        unary_factors = self.unary_factors

        for it in xrange(maxIter):
            for i in nodes:
                i_adj = [x for x in nodes[i]]
                facts = [factors[y][i] for y in i_adj]
                prod = np.prod(facts, axis=0)
                for j in nodes[i]:
                    nodes[i][j] = np.multiply(unary_factors[i], prod)
                    nodes[i][j] = np.divide(nodes[i][j], np.maximum(factors[j][i], [1e-10, 1e-10]))

                    # normalizing
                    nodes[i][j] = np.divide(nodes[i][j], np.sum(nodes[i][j]))

            for j in factors:
                for i in factors[j]:
                    j_adj = [x for x in factors[j] if x != i]
                    node_0 = [nodes[y][j][0] for y in j_adj]
                    node_1 = [nodes[y][j][1] for y in j_adj]
                    diff = np.array(node_0) - np.array(node_1)
                    prod = np.prod(diff)

                    factors[j][i][0] = np.maximum((1.0 + prod)/2, 1e-10)
                    factors[j][i][1] = np.maximum((1.0 - prod)/2, 1e-10)

                    # normalizing
                    factors[j][i] = np.divide(factors[j][i], np.sum(factors[j][i]))
                

            # get marginal probabilities after each iteration
            probs = []
            for i in nodes:
                i_adj = [x for x in nodes[i]]
                facts_0 = [factors[y][i][0] for y in i_adj]
                facts_1 = [factors[y][i][1] for y in i_adj]
                prod_0 = np.maximum(np.prod(facts_0), 1e-10)
                prod_1 = np.maximum(np.prod(facts_1), 1e-10)

                prob = prod_1 / (prod_0 + prod_1)
                probs.append(prob)

            self.mar_props[it] = np.array(probs)
            

        # get beliefs at the end, after all iterations
        probs = []
        for i in nodes:
            i_adj = [x for x in nodes[i]]
            facts_0 = [factors[y][i][0] for y in i_adj]
            facts_1 = [factors[y][i][1] for y in i_adj]
            prod_0 = np.maximum(np.prod(facts_0), 1e-10)
            prod_1 = np.maximum(np.prod(facts_1), 1e-10)

            prob = prod_1 / (prod_0 + prod_1)
            probs.append(prob)

        self.beliefs = np.array(probs)
        ids = self.beliefs < 1e-9
        self.beliefs[ids] = 0

# simulate output with noise
def simulate_output(x, eps):
    out = deepcopy(x)
    ids = np.random.uniform(0, 1, x.shape) < eps
    out[ids] = 1 - out[ids]
    return out

# get realized codeword from probability
def estimate_codeWord(probs):
    est_codeWord = np.zeros((len(probs), 1))
    ids = probs > 0.5
    est_codeWord[ids] = 1
    return est_codeWord

##############################################
##############################################
##############################################

# Simple test case for an invalid codeword
H = np.array([[1,1,0,1,1,1,0,1], [1,0,1,1,1,1,1,0], [1,1,1,1,0,0,1,1], [0,1,1,0,1,1,1,1]])
invalid_codeWord = np.array([1,1,0,1,1,0,1,0])
eps = 0.05
noisy_output = simulate_output(invalid_codeWord, eps)
F = ldpc_graph(H, noisy_output, eps)
F.create_ldpc_graph()

# loopy belief prop

F.loopy_belief_prop(maxIter)
beliefs = F.get_beliefs()

# plot the estimated posterior probability that each codeword bit
# equals one.
nBits = H.shape[1]
xVals = np.linspace(1, nBits, num=nBits)
yVals = beliefs

plt.plot(xVals, yVals, linewidth=2.0)
plt.xlabel('Bit Number')
plt.ylabel('Prob that Each Bit Is One')
plt.title('Estimated Posterior Probability That Each Bit Is One')
plt.show()

