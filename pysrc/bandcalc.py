import tblg_kpy
import numpy as np
import matplotlib.pyplot as plt
import math

kp = tblg_kpy.Kp_tblg_construct()

theta = 1.15

kp.setTwist(theta)
kp.loadFiles("full_relax_kp_01-06-2019.dat")
kp.prepare()

K_pt = kp.getK()
M_pt = kp.getM()
G_pt = kp.getGamma()

num_sections = 3
nk = 20
k_endpoints = np.array([K_pt,G_pt,M_pt,K_pt])

val_list = np.ndarray(shape=(num_sections*nk,kp.getSize()), dtype=float, order='F')
val_list2 = np.ndarray(shape=(num_sections*nk,kp.getSize()), dtype=float, order='F')

k_idx = 0
for k_sec in range(num_sections):

	k_start = k_endpoints[k_sec]
	k_end = k_endpoints[k_sec+1]

	for k0 in range(nk):

		kfac = k0/nk

		k_here = (1-kfac)*k_start + (kfac)*k_end

		H = kp.getH( k_here)
		vals = np.linalg.eigvalsh(H)
		val_list[k_idx,:] = vals

		H2 = kp.getH(-1.0*k_here)
		vals2 = np.linalg.eigvalsh(H2)
		val_list2[k_idx,:] = vals2

		k_idx = k_idx + 1
		print(k_idx/(nk*num_sections))

plt.plot(val_list2,'r')
plt.plot(val_list,'k')
plt.axis([0, num_sections*nk-1, -0.1, 0.1])

plt.show()
