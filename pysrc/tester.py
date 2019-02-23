import tblg_kpy
import numpy as np
import matplotlib.pyplot as plt
import math

kp = tblg_kpy.Kp_tblg_construct()

theta = 1

kp.setTwist(theta)
kp.loadFiles("full_relax_kp_01-06-2019.dat")
kp.prepare()

rot_theta = theta*math.pi/180.0
lattice_a = 1.42*math.sqrt(3.0)

KD = 4*math.pi/(3.0*lattice_a)
KTH = 2*KD*math.sin(rot_theta/2.0)
HEX_BLEN = KTH*math.sqrt(3.0)

hex_b1 = HEX_BLEN*np.array([math.sqrt(3.0)/2.0,-1.0/2.0])
hex_b2 = HEX_BLEN*np.array([0.0,1.0])
hex_shift = (-hex_b1 + hex_b2)/3.0

kk1a = (hex_b1 - hex_b2)/3.0
kk1b = (hex_b1 + 2.0*hex_b2)/3.0
kk1c = (-2.0*hex_b1 - hex_b2)/3.0
kk3 = 0.0*hex_b1
kk4 = hex_b1/2.0

num_sections = 3
nk = 20
k_endpoints = np.array([kk1a,kk3,kk4,kk1a])

val_list = np.ndarray(shape=(num_sections*nk,kp.getSize()), dtype=float, order='F')

k_idx = 0
for k_sec in range(num_sections):

	k_start = k_endpoints[k_sec]
	k_end = k_endpoints[k_sec+1]

	for k0 in range(nk):

		kfac = k0/nk

		k_here = (1-kfac)*k_start + (kfac)*k_end
		H = kp.getH(k_here)
		#print(H_gam.shape)

		vals = np.linalg.eigvalsh(H)
		val_list[k_idx,:] = vals
		k_idx = k_idx + 1
		print(k_idx/(nk*num_sections))

plt.plot(val_list)
plt.show()
