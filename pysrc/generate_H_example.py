import tblg_kpy
import numpy as np
import matplotlib.pyplot as plt
import math

kp = tblg_kpy.Kp_tblg_construct()

theta = 4

kp.setTwist(theta)
kp.loadFiles("full_relax_kp_01-06-2019.dat")
kp.prepare()

K_pt = kp.getK()
M_pt = kp.getM()
G_pt = kp.getGamma()

H = kp.getH(K_pt)
dxH = kp.getGradH(K_pt,0)
dyH = kp.getGradH(K_pt,1)

dk = 1e-8;

K_dx = K_pt.copy();
K_dx[0] = K_dx[0] + dk;

K_dy = K_pt.copy();
K_dy[1] = K_dy[1] + dk;

H_K_dx = kp.getH(K_dx)
H_K_dy = kp.getH(K_dy)

dxH_numeric = (H_K_dx - H)/dk
dyH_numeric = (H_K_dy - H)/dk

print("K_pt: ",K_pt)
print("K_dx: ",K_dx)
print("K_dy: ",K_dy)

#print(dxH_numeric[0:1][0:1])
#print(dxH[0:1][0:1])

print("Max dxH: ",np.max(np.abs(dxH)) )
print("Max dxH_numeric: ", np.max(np.abs(dxH_numeric)) )
print("Max dxH test: ",np.max(np.abs(dxH - dxH_numeric)))
print("Max dyH test: ",np.max(np.abs(dyH - dyH_numeric)))


