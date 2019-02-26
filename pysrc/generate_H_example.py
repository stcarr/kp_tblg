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

H = kp.getH(K_pt)
dxH = kp.getGradH(K_pt,0)
dyH = kp.getGradH(K_pt,1)
