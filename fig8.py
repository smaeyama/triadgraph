#!/usr/bin/env python
# coding: utf-8

# # 1. Replot time evolution of kinetic energy, Fig 5 (a) in Sasaki, PPCF (2020).

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rcParams['figure.figsize'] = (2.3,4.2) # 幅,高さ(inch)
plt.rcParams["figure.dpi"] = 150       # dpi(dots per inch)

plt.rcParams['font.size'] = 7 #フォントの大きさ

plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')

plt.rcParams["legend.edgecolor"] = "k"

plt.rcParams["axes.edgecolor"] = "k"
plt.rcParams['axes.linewidth'] = 0.8

#plt.tight_layout()




####
data_A=np.loadtxt("./sample_data/J_evo_wo_sym_A.dat")
data_B=np.loadtxt("./sample_data/J_evo_wo_sym_B.dat")
data_C=np.loadtxt("./sample_data/J_evo_wo_sym_C.dat")
data_D=np.loadtxt("./sample_data/J_evo_wo_sym_D.dat")

nodename=["A","B","C","D"]
#nodename=["Background_Flow_Shear","Zonal_Flow","Kelvin-Helmholtz_Mode","Intermittent_Spiral_Structure"]
#nodename=["B","Z","K","S"]
a_kpq = np.array([data_A.reshape(4,4,1001),data_B.reshape(4,4,1001),data_C.reshape(4,4,1001),data_D.reshape(4,4,1001)])
a_kpq=a_kpq[:,:,:,120:221] # *** Time slice ***
from triadgraph import symmetrize_triadtransfer, directional_triadtransfer,                          triadgraph_symmetric_kpq, triadgraph_symmetric_all,                          triadgraph_directional_kpq, triadgraph_directional_all
j_kpq = symmetrize_triadtransfer(a_kpq)
d_kpq=directional_triadtransfer(j_kpq)
# Normalization for plot
#a_kpq_max = np.max(abs(a_kpq))
j_kpq_max = np.max(abs(j_kpq))
#d_kpq_max = np.max(abs(d_kpq))
a_kpq = a_kpq / j_kpq_max
j_kpq = j_kpq / j_kpq_max
d_kpq = d_kpq / j_kpq_max
wa_kpq = np.transpose(a_kpq,axes=(0,2,1,3))
####




fig = plt.figure()
ax = fig.add_subplot(211)
surf = ax.pcolormesh(np.arange(4),np.arange(4),np.sum(np.average(d_kpq[:,:,:,50:60],axis=3),axis=1),
                     cmap=cm.RdBu,shading="auto",vmax=0.7,vmin=-0.7)
ax.set_xlabel("Mode k")
ax.set_ylabel("Mode q")
ax.set_xticks(list(np.arange(4)))
ax.set_xticklabels(nodename)
ax.set_yticks(list(np.arange(4)))
ax.set_yticklabels(nodename)
ax.text(-0.4,3.15,'(a)')
ax.set_aspect('equal', adjustable='box')
cbar = fig.colorbar(surf,aspect=5)
cbar.set_label(r"$\sum_p D_{k \leftarrow q}^p$")

ax = fig.add_subplot(212)
surf = ax.pcolormesh(np.arange(4),np.arange(4),np.sum(np.average(wa_kpq[:,:,:,50:60],axis=3),axis=1),
                     cmap=cm.RdBu,shading="auto",vmax=0.7,vmin=-0.7)
ax.set_xlabel("Mode k")
ax.set_ylabel("Mode q")
ax.set_xticks(list(np.arange(4)))
ax.set_xticklabels(nodename)
ax.set_yticks(list(np.arange(4)))
ax.set_yticklabels(nodename)
ax.text(-0.4,3.15,'(b)')
ax.set_aspect('equal', adjustable='box')
cbar = fig.colorbar(surf,aspect=5)
cbar.set_label(r"$\sum_p A_k^{q,p}$")


# In[ ]:





# In[ ]:




fig.savefig("fig8.pdf", bbox_inches="tight")
#plt.show()

plt.close(fig)

