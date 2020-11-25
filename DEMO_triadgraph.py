#!/usr/bin/env python
# coding: utf-8

# # 1. Replot time evolution of kinetic energy, Fig 5 (a) in Sasaki, PPCF (2020).

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rcParams["font.size"]=14 # Font size for matplotlib
#plt.tight_layout()

energy=np.loadtxt("./sample_data/Energy_evo.dat")

nodename=["A","B","C","D"]
#nodename=["Background_Flow_Shear","Zonal_Flow","Kelvin-Helmholtz_Mode","Intermittent_Spiral_Structure"]
#nodename=["B","Z","K","S"]

### Plot Energy
fig=plt.figure(figsize=[6,2])
ax = fig.add_subplot(111)
for k in range(4):
    ax.plot(energy[k,:],label="$E_{" + nodename[k] + "}$")
ax.set_xlabel("Time t")
ax.set_ylabel("Energy $E_k$")
ax.set_ylim(0,None)
ax.set_xlim(0,None)
ax.legend(bbox_to_anchor=(1.01,0.5), loc="center left", borderaxespad=0, labelspacing=0)
plt.show()


# # 2. Replot time evolution of transfer, Fig 6 (e)-(h) in Sasaki, PPCF (2020).

# In[ ]:


data_A=np.loadtxt("./sample_data/J_evo_wo_sym_A.dat")
data_B=np.loadtxt("./sample_data/J_evo_wo_sym_B.dat")
data_C=np.loadtxt("./sample_data/J_evo_wo_sym_C.dat")
data_D=np.loadtxt("./sample_data/J_evo_wo_sym_D.dat")

vmin=-0.02
vmax=0.02

fig = plt.figure(figsize=(8,8))

ax = fig.add_subplot(411)
surf = ax.pcolormesh(data_A,cmap=cm.jet,vmin=vmin,vmax=vmax)
ax.set_xlabel("Time t")
ax.set_ylabel("Coupling ID")
cbar = fig.colorbar(surf,aspect=5)
cbar.set_label(nodename[0])

ax = fig.add_subplot(412)
surf = ax.pcolormesh(data_B,cmap=cm.jet,vmin=vmin,vmax=vmax)
ax.set_xlabel("Time t")
ax.set_ylabel("Coupling ID")
cbar = fig.colorbar(surf,aspect=5)
cbar.set_label(nodename[1])

ax = fig.add_subplot(413)
surf = ax.pcolormesh(data_C,cmap=cm.jet,vmin=vmin,vmax=vmax)
ax.set_xlabel("Time t")
ax.set_ylabel("Coupling ID")
cbar = fig.colorbar(surf,aspect=5)
cbar.set_label(nodename[2])

ax = fig.add_subplot(414)
surf = ax.pcolormesh(data_D,cmap=cm.jet,vmin=vmin,vmax=vmax)
ax.set_xlabel("Time t")
ax.set_ylabel("Coupling ID")
cbar = fig.colorbar(surf,aspect=5)
cbar.set_label(nodename[3])

plt.show()
#fig.savefig(file)
plt.close(fig)


# # 3. Rearrange data of energy transfer, and construct non-symmetrized triad transfer $A_k^{p,q}$
# Anti-symmetry $A_k^{p,q} = - A_p^{k,q}$

# In[ ]:


# Rearrange data
# 0: A - Background_Flow_Shear
# 1: B - Zonal_Flow
# 2: C - Kelvin-Helmholtz_Mode
# 3: D - Intermittent_Spiral_Structure
# Assymmetric triad transfer A_k^pq

a_kpq = np.array([data_A.reshape(4,4,1001),data_B.reshape(4,4,1001),data_C.reshape(4,4,1001),data_D.reshape(4,4,1001)])


energy=energy[:,120:221]   # *** Time slice ***
a_kpq=a_kpq[:,:,:,120:221] # *** Time slice ***



# for k in range(4):
#     for p in range(4):
#         for q in range(4):
#             if (k!=p):
#                 print(k,p,q,a_kpq[k,p,q,30],a_kpq[p,k,q,30],(a_kpq[k,p,q,30]+a_kpq[p,k,q,30])/max(abs(a_kpq[k,p,q,30]),abs(a_kpq[p,k,q,30])))

print(a_kpq[0,1,2,30],a_kpq[1,0,2,30], "# J(zeta|alpha,beta) = -J(alpha|zeta,beta) in Eq. (5b) in Sasaki PoP")
print(a_kpq[0,1,1,30],a_kpq[1,0,1,30], "# J(zeta|alpha,beta) = -J(alpha|zeta,beta) in Eq. (5b) in Sasaki PoP")
print(a_kpq[2,3,3,30],a_kpq[3,2,3,30], "# Combination (2,3,3) has a bit large error.")


# # 4. Symmetrization from non-symmetrized triad transfer $A_k^{p,q}$ to the symmetrized triad transfer $S_k^{p,q}$
# $S_k^{p,q}= \frac{1}{2} (A_k^{p,q} + A_k^{q,p})$  
# Symmetry $S_k^{p,q} = S_k^{q,p}$  
# Detailed balance $S_k^{p,q} + S_p^{q,k} S_q^{k,p} = 0$

# In[ ]:


from triadgraph import symmetrize_triadtransfer, directional_triadtransfer,                          triadgraph_symmetric_kpq, triadgraph_symmetric_all,                          triadgraph_directional_kpq, triadgraph_directional_all

j_kpq = symmetrize_triadtransfer(a_kpq)

print(j_kpq[0,1,2,30],j_kpq[0,2,1,30])
print(j_kpq[0,1,2,30],j_kpq[1,2,0,30],j_kpq[2,0,1,30],j_kpq[0,1,2,30]+j_kpq[1,2,0,30]+j_kpq[2,0,1,30])
print(j_kpq[2,3,3,30],j_kpq[3,2,3,30],j_kpq[3,3,2,30],j_kpq[2,3,3,30]+j_kpq[3,2,3,30]+j_kpq[3,3,2,30])


# # 5. Plot the net transfer $T_k = \sum_p \sum_q S_k^{p,q}$, and replot time evolution of energy transfer, Fig. 7 (b)-(e) in Sasaki, PPCF (2020).

# In[ ]:


t_k = np.sum(np.sum(j_kpq,axis=2),axis=1)

### Plot all T_k
fig=plt.figure(figsize=[6,2])
ax = fig.add_subplot(111)
for k in range(4):
    ax.plot(t_k[k,:],label="$T_{" + nodename[k] + "}$")
ax.plot(np.sum(t_k,axis=0),label="sum")
ax.set_xlabel("Time t")
ax.set_ylabel("Net transfer $T_k$")
ax.set_xlim(0,100)
ax.legend(bbox_to_anchor=(1.01,0.5), loc="center left", borderaxespad=0, labelspacing=0)
plt.show()

### Plot J_k^p,q for each k
fig = plt.figure(figsize=[6,9])
screening = 0.1 * np.max(abs(j_kpq))
for k in range(4):
    ax = fig.add_subplot(410+(k+1))
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
    ax.set_ylabel("$J_{"+ nodename[k] + "}^{p,q}$")
    ax.set_xlim(0,100)
    ax.legend(bbox_to_anchor=(1.01,0.5), loc="center left", borderaxespad=0, labelspacing=0)
plt.show()

plt.close(fig)


# # 6. Directional representation $D_{k \leftarrow q}^p$ constructed from the symmetrized triad transfer $S_k^{p,q}$

# In[ ]:


d_kpq=directional_triadtransfer(j_kpq)
print(j_kpq[0,1,2,10],j_kpq[1,2,0,10],j_kpq[2,0,1,10],j_kpq[0,1,2,10]+j_kpq[1,2,0,10]+j_kpq[2,0,1,10])
print(d_kpq[0,1,2,10],d_kpq[2,1,0,10])
print(d_kpq[1,2,0,10],d_kpq[0,2,1,10])
print(d_kpq[2,0,1,10],d_kpq[1,0,2,10])

print(j_kpq[0,2,2,10],j_kpq[2,2,0,10],j_kpq[2,0,2,10],j_kpq[0,2,2,10]+j_kpq[2,2,0,10]+j_kpq[2,0,2,10])
print(d_kpq[0,2,2,10],d_kpq[2,2,0,10])
print(d_kpq[2,2,0,10],d_kpq[0,2,2,10])
print(d_kpq[2,0,2,10],d_kpq[2,0,2,10])

print(np.sum(np.sum(a_kpq[:,:,:,10],axis=2),axis=1))
print(np.sum(np.sum(j_kpq[:,:,:,10],axis=2),axis=1))
print(np.sum(np.sum(d_kpq[:,:,:,10],axis=2),axis=1))


# # 7. Network visualization
# Since the relative amplitude of energy transfer is expressed by the width of arrows, normalization is recommended before the plot.

# In[ ]:


# Normalization for plot
#a_kpq_max = np.max(abs(a_kpq))
j_kpq_max = np.max(abs(j_kpq))
#d_kpq_max = np.max(abs(d_kpq))
a_kpq = a_kpq / j_kpq_max
j_kpq = j_kpq / j_kpq_max
d_kpq = d_kpq / j_kpq_max
energy_max = np.max(abs(energy))
energy = energy / energy_max


# In[ ]:


triadgraph_symmetric_kpq(j_kpq[:,:,:,72],1,2,3,title="(a) symmetric triad transfer among three",pwidth=10,nodename=nodename)
triadgraph_directional_kpq(d_kpq[:,:,:,72],1,2,3,title="(b) directional representation of (a)",pwidth=10,nodename=nodename)

triadgraph_symmetric_kpq(j_kpq[:,:,:,72],2,3,3,title="(c) symmetric triad transfer between two",pwidth=5,nodename=nodename)
triadgraph_directional_kpq(d_kpq[:,:,:,72],2,3,3,title="(d) directional representation of (c)",pwidth=5,nodename=nodename)


# ## 7.1 Visualization of symmetrized triad transfer $S_k^{p,q}$

# In[ ]:


for i in range(j_kpq.shape[3]):
    triadgraph_symmetric_all(j_kpq[:,:,:,i],output="png/symmetric_jkpq_t{0:08d}.png".format(i),title="t={:4d}".format(i),nodename=nodename,energy=energy[:,i])
    triadgraph_symmetric_all(j_kpq[:,:,:,i],output="dot/symmetric_jkpq_t{0:08d}.dot".format(i),title="t={:4d}".format(i),nodename=nodename,energy=energy[:,i])
for i in range(10):
    triadgraph_symmetric_all(np.average(j_kpq[:,:,:,i*10:(i+1)*10],axis=3),title="t={:4d}-{:4d}".format(i*10,(i+1)*10),nodename=nodename,energy=np.average(energy[:,i*10:(i+1)*10],axis=1))


# ## 7.2 Visualization of directional representation $D_{k \leftarrow q}^p$
# We propose the directional representation which looks qualitatively similar to the symmetrized triad transfer.

# In[ ]:


for i in range(d_kpq.shape[3]):
    triadgraph_directional_all(d_kpq[:,:,:,i],output="png/directional_dkpq_t{0:08d}.png".format(i),title="t={:4d}".format(i),nodename=nodename,energy=energy[:,i])
    triadgraph_directional_all(d_kpq[:,:,:,i],output="dot/directional_dkpq_t{0:08d}.dot".format(i),title="t={:4d}".format(i),nodename=nodename,energy=energy[:,i])
for i in range(10):
    triadgraph_directional_all(np.average(d_kpq[:,:,:,i*10:(i+1)*10],axis=3),title="t={:4d}-{:4d}".format(i*10,(i+1)*10),nodename=nodename,energy=np.average(energy[:,i*10:(i+1)*10],axis=1))


# ## 7.3 Visualization of non-symmetrized triad transfer $A_k^{p,q}$
# One may observe fictitious interactions and qualitatively different interpretation from the symmetrized transfer.

# In[ ]:


wa_kpq = np.transpose(a_kpq,axes=(0,2,1,3))
for i in range(wa_kpq.shape[3]):
    triadgraph_directional_all(wa_kpq[:,:,:,i],output="png/asymmetric_akpq_t{0:08d}.png".format(i),title="t={:4d}".format(i),nodename=nodename,energy=energy[:,i])
    triadgraph_directional_all(wa_kpq[:,:,:,i],output="dot/asymmetric_akpq_t{0:08d}.dot".format(i),title="t={:4d}".format(i),nodename=nodename,energy=energy[:,i])
for i in range(10):
    triadgraph_directional_all(np.average(wa_kpq[:,:,:,i*10:(i+1)*10],axis=3),title="t={:4d}-{:4d}".format(i*10,(i+1)*10),nodename=nodename,energy=np.average(energy[:,i*10:(i+1)*10],axis=1))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:



fig=plt.figure(figsize=[6,11])
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
plt.show()

plt.close(fig)


# In[ ]:


# comparison of visivility
triadgraph_symmetric_all(np.average(j_kpq[:,:,:,50:60],axis=3),title="t={:4d}-{:4d}".format(50,60),nodename=nodename,screening=0,energy=np.average(energy[:,i*10:(i+1)*10],axis=1))
triadgraph_directional_all(np.average(d_kpq[:,:,:,50:60],axis=3),title="t={:4d}-{:4d}".format(50,60),nodename=nodename,screening=0,energy=np.average(energy[:,i*10:(i+1)*10],axis=1))


# In[ ]:


fig = plt.figure()
ax = fig.add_subplot(111)
surf = ax.pcolormesh(np.arange(4),np.arange(4),np.sum(np.average(d_kpq[:,:,:,50:60],axis=3),axis=1),
                     cmap=cm.RdBu,shading="auto",vmax=0.7,vmin=-0.7)
ax.set_xlabel("Mode k")
ax.set_ylabel("Mode q")
ax.set_xticks(list(np.arange(4)))
ax.set_xticklabels(nodename)
ax.set_yticks(list(np.arange(4)))
ax.set_yticklabels(nodename)
cbar = fig.colorbar(surf,aspect=5)
cbar.set_label(r"$\sum_p D_{k \leftarrow q}^p$")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
surf = ax.pcolormesh(np.arange(4),np.arange(4),np.sum(np.average(wa_kpq[:,:,:,50:60],axis=3),axis=1),
                     cmap=cm.RdBu,shading="auto",vmax=0.7,vmin=-0.7)
ax.set_xlabel("Mode k")
ax.set_ylabel("Mode q")
ax.set_xticks(list(np.arange(4)))
ax.set_xticklabels(nodename)
ax.set_yticks(list(np.arange(4)))
ax.set_yticklabels(nodename)
cbar = fig.colorbar(surf,aspect=5)
cbar.set_label(r"$\sum_p A_k^{q,p}$")
plt.show()


# In[ ]:





# In[ ]:




