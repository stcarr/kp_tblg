import tblg_kpy
import numpy as np
import matplotlib.pyplot as plt
import math

kp = tblg_kpy.Kp_tblg_construct()

theta = 1

kp.setTwist(theta)
kp.loadFiles("full_relax_kp_01-06-2019.dat")
kp.prepare()

G_pt = kp.getGamma()

H = kp.getH(G_pt)

index_to_G = kp.getIndexToG()
G_to_index = kp.getGToIndex()

# the search size we use for finding connections
hex_M = 100;

# G_to_index is of size 2*hex_M+1  x 2*hex_M+1
# [G1,G2] = [n1,n2] is G_to_index[n1+100,n2+100]

# number of unique [G1,G2] for one layer:
num_hex = index_to_G.size/2

# hex_shift = (-G1 + G2)/3
# momentum is given by:
# Layer 1:
# G1*n1 + G2*n2 + hex_shift
# Layer 2:
# -(G1*n1 + G2*n2 + hex_shift)

# if we are at K given by Layer 1: n1,n2
# the three nearest neighbors on Layer 2 given by:
# -[n1, n2+1]
# -[n1-1,n2+1]
# -[n1-1,n2]

# indexing for H goes as follows:
# Index in H    : n,            Layer # and orbital (A/B)
# 1             : 0,            L1A
# 2             : 0,            L1B
# 3             : 1,            L1A
# 4             : 1,            L1B
# ...
# 2*num_hex-1   : num_hex-1,    L1A
# 2*num_hex     : num_hex-1,    L1B
# 2*num_hex + 1 : 0,            L2A
# 2*num_hex + 2 : 0,            L2B
# 2*num_hex + 3 : 1,            L2A
# 2*num_hex + 4 : 1,            L2B
# ... etc


# lets find all nearest-neighbor couplings to an element of our matrix
# EVERYTHING IS BASE 0 !
L1_k_index = G_to_index[0+hex_M,0+hex_M]
L1_orb = 0 # A or B orbital on layer 1 [A = 0, B = 1]
L2_orb = 1 # A or B orbital on layer 2

# Element of the Hamiltonian
L1_H_index = 2*L1_k_index + L1_orb # orbital of the k_index'th k point

L1_n1 = index_to_G[L1_k_index][0]
L1_n2 = index_to_G[L1_k_index][1]

# The three nearest-neighbor indices
# minus sign as layer 1 does:
# k_[n1,n2] = +(G1*n1 + G2*n2 + hex_shift)
# whike layer 2 does:
# k_[n1,n2] = -(G1*n1 + G2*n2 + hex_shift)
L2_n1_1 = -(L1_n1+0)
L2_n2_1 = -(L1_n2+1)

L2_n1_2 = -(L1_n1-1)
L2_n2_2 = -(L1_n2+1)

L2_n1_3 = -(L1_n1-1)
L2_n2_3 = -(L1_n2+0)

# get the index of the k point
L2_k_index_1 = G_to_index[L2_n1_1+hex_M ,L2_n2_1+hex_M]
L2_k_index_2 = G_to_index[L2_n1_2+hex_M ,L2_n2_2+hex_M]
L2_k_index_3 = G_to_index[L2_n1_3+hex_M ,L2_n2_3+hex_M]

# get the index in the Hamiltonian
L2_H_index_1 = 2*num_hex + 2*L2_k_index_1 + L2_orb
L2_H_index_2 = 2*num_hex + 2*L2_k_index_2 + L2_orb
L2_H_index_3 = 2*num_hex + 2*L2_k_index_3 + L2_orb

# Check the output
print("N.N. Coupling 1: ",H[L1_H_index,L2_H_index_1])
print("N.N. Coupling 2: ",H[L1_H_index,L2_H_index_2])
print("N.N. Coupling 3: ",H[L1_H_index,L2_H_index_3])

# Note that we use the convention that the k-pts for Layer 2
# are equal to the negative of the kpts of Layer 1.
# for example: if n = 0 gives [G1,G2] = [5,4]
# then for layer 1: This is at momenta G1*5 + G2*4
# but for layer  2: This is at momenta G1*-5 + G2*-4
# So if n1 gave G1,G2] = [5,0] and n2 gave [G1,G2] = [-5,0]
# we would expect 2*n1 to couple to 2*num_hex+2*n2 in the Hamiltonian
