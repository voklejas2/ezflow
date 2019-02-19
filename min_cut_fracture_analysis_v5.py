import fracture_analysis_vo_mod_v2 as fa
import networkx as nx
from networkx.algorithms.flow import shortest_augmenting_path
from networkx.algorithms.flow import edmonds_karp
import numpy as np
import pickle
import matplotlib.pyplot as plt


#bond_file = './data170511/topology_symmetric_592.txt'
#coord_file = './data170511/coords_symmetric_592.txt'
#bond_file = 'topology_587_delete.txt'
#coord_file = 'coords_delete_587.txt'
bond_file = 'topology_burg.txt'
coord_file = 'coords_burg.txt'


direction = 'x_coord'
attribute = 'rupture_force'
algorithm = shortest_augmenting_path

sim = fa.FractureAnalysis(bond_file, coord_file)
sim.find_center(direction)
print sim.min_val, sim.max_val

cut_values = []
cut_edges = []
cut_freq = []
cut_height = []
cut_edges_all = []
for delta in np.arange(0.5,25.0,0.1):
    print 'delta = %2.3f' % delta
    sim.update_bond_coord_lists(direction, delta)
    sim.generate_directed_graph()
    print 'running cut algorithm...'
    cut_value, edge_cut_list = sim.min_cut(direction, attribute, delta, algorithm)
    print 'min cut completed...'
    cut_values.append([delta, cut_value])
    height = sim.find_fracture_height(edge_cut_list, direction)
    print delta, height
    cut_height.append([delta, height])
    for edge in edge_cut_list:
        if edge in cut_edges:
            ind = cut_edges.index(edge)
            cut_freq[ind] += 1
        else:
            cut_edges.append(edge)
            cut_freq.append(1)
    cut_edges_all.append([delta, edge_cut_list])
    print 'data lists appended...'
    print ''
#    print nx.get_edge_attributes(sim.g,'rupture_force')
#    print cut_value
    sim.reset_graph()

#pickle.dump(cut_values, open('./data170511/sym_592_cutvalue_ap-z.p','wb'))
#pickle.dump(cut_height, open('./data170511/sym_592_fractureheight_ap-z.p','wb'))
#pickle.dump(cut_freq, open('./data170511/sym_592_edgefreq_ap-z.p','wb'))
#pickle.dump(cut_edges_all, open('./data170511/sym_592_cutedges_ap-z.p','wb'))


pickle.dump(cut_values, open('burg_cutvalue_ap-x.p','wb'))
pickle.dump(cut_height, open('burg_fractureheight_ap-x.p','wb'))
pickle.dump(cut_freq, open('burg_edgefreq_ap-x.p','wb'))
pickle.dump(cut_edges_all, open('burg_cutedges_ap-x.p','wb'))


#pickle.dump(cut_values, open('center_void_quartz_cutvalue_ap-y.p','wb'))
#pickle.dump(cut_height, open('center_void_quartz_fractureheight_ap-y.p','wb'))
#pickle.dump(cut_freq, open('center_void_quartz_edgefreq_ap-y.p','wb'))
#pickle.dump(cut_edges_all, open('center_void_quartz_cutedges_ap-y.p','wb'))

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.plot(cut_values[:,0], cut_values[:,1], '-bo')
plt.xlabel('Delta (Angstrom)', fontsize=12)
plt.ylabel('Edge Weight Threshold (Stress, GPa)', fontsize=12)  

cut_height = np.array(cut_height)
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.plot(cut_height[:,0], cut_height[:,1], '-bo')
plt.xlabel('Delta (Angstrom)', fontsize=12)
plt.ylabel('Fracture Height (Angstrom)', fontsize=12)  
plt.show()
print "stop\n"
print "print is not a function?"