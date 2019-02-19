import numpy as np
import matplotlib.pyplot as plt
import pickle





#-----------------------------------------------------
#    Plot rupture force and fracture height vs. Delta
#-----------------------------------------------------

cut_values = pickle.load(open('center_void_quartz_cutvalue_ek-z.p','rb'))
cut_height = pickle.load(open('center_void_quartz_fractureheight_ek-z.p','rb'))
cut_freq = pickle.load(open('center_void_quartz_edgefreq_ek-z.p','rb'))
cut_edges_all = pickle.load(open('center_void_quartz_cutedges_ek-z.p','rb'))

fig1 = plt.figure(dpi=200)

# Rupture force
cut_values = np.array(cut_values)
ax = fig1.add_subplot(121)
ax.plot(2.0*cut_values[:,0], cut_values[:,1], '-b.')
ax.set_xlabel(r'Delta ($\AA$)', fontsize=13)
ax.set_ylabel('Total Rupture Force (nN)', fontsize=13)

# Rupture force - inset plot
cut_values = np.array(cut_values)
left, bottom, width, height = [0.30,0.22,0.16,0.2]
ax = fig1.add_axes([left, bottom, width, height])
#ax = fig1.add_subplot(122)
#ax.plot(2.0*cut_values[:,0][60:150], cut_values[:,1][60:150], '-b.')
ax.plot(2.0*cut_values[:,0][0:38], cut_values[:,1][0:38], '-b.')
ax.tick_params(axis='both', which='major', labelsize=5)
ax.set_xlim([1.0,3.8])
ax.set_ylim([0,450])
ax.set_xlabel(r'Delta ($\AA$)', fontsize=7)
ax.set_ylabel('Total Rupture Force (nN)', fontsize=7)


# Fracture height
cut_height = np.array(cut_height)
ax = fig1.add_subplot(122)
ax.plot(2.0*cut_height[:,0], cut_height[:,1], '-g.')
ax.set_xlabel(r'Delta ($\AA$)', fontsize=13)
ax.set_ylabel(r'Fracture Height ($\AA$)', fontsize=13)
# Rupture force - line
lin_x = np.linspace(0,2.0*cut_height[:,0].max(),100)
ax.plot(lin_x,lin_x,'--r')

# Fracture height - inset plot
cut_height = np.array(cut_height)
left, bottom, width, height = [0.78,0.22,0.16,0.2]
ax = fig1.add_axes([left, bottom, width, height])
#ax = fig1.add_subplot(124)
ax.plot(2.0*cut_height[:,0], cut_height[:,1], '-g.')
lin_x = np.linspace(1.3,2.0*cut_height[:,0].max(),100)
ax.plot(lin_x,lin_x,'--r')
ax.tick_params(axis='both', which='major', labelsize=5)
ax.set_xlim([1.0,3.8])
ax.set_ylim([1.2,2.5])
ax.set_xlabel(r'Delta ($\AA$)', fontsize=7)
ax.set_ylabel(r'Fracture Height ($\AA$)', fontsize=7)
#plt.tight_layout()
plt.show()


#-----------------------------------------------------
#    Plot frequency of edges in min cut
#-----------------------------------------------------
fig2 = plt.figure()
cut_freq = np.array(cut_freq)
ax = fig2.add_subplot(111)
ax.hist(cut_freq,100)
ax.set_xlabel('Cut Edges', fontsize=15)
ax.set_ylabel('Freq', fontsize=15)
#plt.tight_layout()


