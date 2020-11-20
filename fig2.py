#!/usr/bin/env python
# coding: utf-8

# # Fig 2 in Maeyama, NJP (2020).

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = (3.14, 3.14)
plt.rcParams['font.family'] ='sans-serif'#使用するフォント
plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['font.size'] = 8 #フォントの大きさ
plt.rcParams["legend.edgecolor"] = "k"


# In[ ]:


energy=np.loadtxt("./sample_data/Energy_evo.dat")

nodename=["A","B","C","D"]
#nodename=["Background_Flow_Shear","Zonal_Flow","Kelvin-Helmholtz_Mode","Intermittent_Spiral_Structure"]
#nodename=["B","Z","K","S"]

data_A=np.loadtxt("./sample_data/J_evo_wo_sym_A.dat")
data_B=np.loadtxt("./sample_data/J_evo_wo_sym_B.dat")
data_C=np.loadtxt("./sample_data/J_evo_wo_sym_C.dat")
data_D=np.loadtxt("./sample_data/J_evo_wo_sym_D.dat")

# Rearrange data
# 0: A - Background_Flow_Shear
# 1: B - Zonal_Flow
# 2: C - Kelvin-Helmholtz_Mode
# 3: D - Intermittent_Spiral_Structure
# Assymmetric triad transfer A_k^pq
a_kpq = np.array([data_A.reshape(4,4,1001),data_B.reshape(4,4,1001),data_C.reshape(4,4,1001),data_D.reshape(4,4,1001)])


energy=energy[:,120:221]   # *** Time slice ***
a_kpq=a_kpq[:,:,:,120:221] # *** Time slice ***




from triadgraph import symmetrize_triadtransfer, directional_triadtransfer,                          triadgraph_symmetric_kpq, triadgraph_symmetric_all,                          triadgraph_directional_kpq, triadgraph_directional_all
# Symmetrized triad transfer S_k^pq
j_kpq = symmetrize_triadtransfer(a_kpq)


# Net transfer T_k
t_k = np.sum(np.sum(j_kpq,axis=2),axis=1)


# In[ ]:



fig=plt.figure(figsize=[2,6])

ax = fig.add_subplot(611)
for k in range(4):
    ax.plot(energy[k,:],label="$E_{" + nodename[k] + "}$")
ax.set_xlabel("Time t")
ax.set_ylabel("Energy $E_k$")
ax.set_ylim(0,None)
ax.set_xlim(0,100)
plt.xticks(color="None")
ax.text(2,0.67,'(a)')
ax.legend(bbox_to_anchor=(1.03,0.5), loc="center left", borderaxespad=0, labelspacing=0)

ax = fig.add_subplot(612)
for k in range(4):
    ax.plot(t_k[k,:],label="$T_{" + nodename[k] + "}$")
ax.plot(np.sum(t_k,axis=0),label="sum")
ax.set_xlabel("Time t")
ax.set_ylabel("Net transfer $T_k$")
ax.set_xlim(0,100)
plt.xticks(color="None")
ax.text(2,0.067,'(b)')
ax.legend(bbox_to_anchor=(1.03,0.5), loc="center left", borderaxespad=0, labelspacing=0)

screening = 0.1 * np.max(abs(j_kpq))
k=0
ax = fig.add_subplot(612+(k+1))
for p in range(4):
    for q in range(p,4):
        if (p!=q): # J_k^kq = J_k^qk,  J_k^pk = J_k^kp
            if np.max(abs(j_kpq[k,p,q,:])*2) > screening: # Screening for visibility
                #ax.plot(j_kpq[k,p,q,:]*2,label="J({}|{},{})*2".format(nodename[k],nodename[p],nodename[q]))
                ax.plot(j_kpq[k,p,q,:]*2,label="2$J_{" + nodename[k] + "}^{" + nodename[p] + "," + nodename[q] + "}$")
        else: # J_k^pq
            if np.max(abs(j_kpq[k,p,q,:])) > screening: # Screening for visibility
                #ax.plot(j_kpq[k,p,q,:],label="J({}|{},{})".format(nodename[k],nodename[p],nodename[q]))
                ax.plot(j_kpq[k,p,q,:],label="$J_{" + nodename[k] + "}^{" + nodename[p] + "," + nodename[q] + "}$")
ax.set_xlabel("Time t")
ax.set_ylabel("Triad $J_{"+ nodename[k] + "}^{p,q}$")
ax.set_xlim(0,100)
plt.xticks(color="None")
ax.text(2,0.013,'(c)')
ax.legend(bbox_to_anchor=(1.03,0.5), loc="center left", borderaxespad=0, labelspacing=0)

k=1
ax = fig.add_subplot(612+(k+1))
for p in range(4):
    for q in range(p,4):
        if (p!=q): # J_k^kq = J_k^qk,  J_k^pk = J_k^kp
            if np.max(abs(j_kpq[k,p,q,:])*2) > screening: # Screening for visibility
                #ax.plot(j_kpq[k,p,q,:]*2,label="J({}|{},{})*2".format(nodename[k],nodename[p],nodename[q]))
                ax.plot(j_kpq[k,p,q,:]*2,label="2$J_{" + nodename[k] + "}^{" + nodename[p] + "," + nodename[q] + "}$")
        else: # J_k^pq
            if np.max(abs(j_kpq[k,p,q,:])) > screening: # Screening for visibility
                #ax.plot(j_kpq[k,p,q,:],label="J({}|{},{})".format(nodename[k],nodename[p],nodename[q]))
                ax.plot(j_kpq[k,p,q,:],label="$J_{" + nodename[k] + "}^{" + nodename[p] + "," + nodename[q] + "}$")
ax.set_xlabel("Time t")
ax.set_ylabel("Triad $J_{"+ nodename[k] + "}^{p,q}$")
ax.set_xlim(0,100)
plt.xticks(color="None")
ax.text(2,0.038,'(d)')
ax.legend(bbox_to_anchor=(1.03,0.5), loc="center left", borderaxespad=0, labelspacing=0)

k=2
ax = fig.add_subplot(612+(k+1))
for p in range(4):
    for q in range(p,4):
        if (p!=q): # J_k^kq = J_k^qk,  J_k^pk = J_k^kp
            if np.max(abs(j_kpq[k,p,q,:])*2) > screening: # Screening for visibility
                #ax.plot(j_kpq[k,p,q,:]*2,label="J({}|{},{})*2".format(nodename[k],nodename[p],nodename[q]))
                ax.plot(j_kpq[k,p,q,:]*2,label="2$J_{" + nodename[k] + "}^{" + nodename[p] + "," + nodename[q] + "}$")
        else: # J_k^pq
            if np.max(abs(j_kpq[k,p,q,:])) > screening: # Screening for visibility
                #ax.plot(j_kpq[k,p,q,:],label="J({}|{},{})".format(nodename[k],nodename[p],nodename[q]))
                ax.plot(j_kpq[k,p,q,:],label="$J_{" + nodename[k] + "}^{" + nodename[p] + "," + nodename[q] + "}$")
ax.set_xlabel("Time t")
ax.set_ylabel("Triad $J_{"+ nodename[k] + "}^{p,q}$")
ax.set_xlim(0,100)
plt.xticks(color="None")
ax.text(2,0.009,'(e)')
ax.legend(bbox_to_anchor=(1.03,0.5), loc="center left", borderaxespad=0, labelspacing=0)

k=3
ax = fig.add_subplot(612+(k+1))
for p in range(4):
    for q in range(p,4):
        if (p!=q): # J_k^kq = J_k^qk,  J_k^pk = J_k^kp
            if np.max(abs(j_kpq[k,p,q,:])*2) > screening: # Screening for visibility
                #ax.plot(j_kpq[k,p,q,:]*2,label="J({}|{},{})*2".format(nodename[k],nodename[p],nodename[q]))
                ax.plot(j_kpq[k,p,q,:]*2,label="2$J_{" + nodename[k] + "}^{" + nodename[p] + "," + nodename[q] + "}$")
        else: # J_k^pq
            if np.max(abs(j_kpq[k,p,q,:])) > screening: # Screening for visibility
                #ax.plot(j_kpq[k,p,q,:],label="J({}|{},{})".format(nodename[k],nodename[p],nodename[q]))
                ax.plot(j_kpq[k,p,q,:],label="$J_{" + nodename[k] + "}^{" + nodename[p] + "," + nodename[q] + "}$")
ax.set_xlabel("Time t")
ax.set_ylabel("Triad $J_{"+ nodename[k] + "}^{p,q}$")
ax.set_xlim(0,100)
plt.xticks(color="None")
ax.text(2,0.022,'(f)')
ax.legend(bbox_to_anchor=(1.03,0.5), loc="center left", borderaxespad=0, labelspacing=0)

plt.xticks(color="k")


fig.savefig("fig2.pdf", bbox_inches="tight")
plt.show()

plt.close(fig)


# In[ ]:




