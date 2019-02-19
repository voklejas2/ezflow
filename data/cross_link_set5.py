import fracture_analysis_vo_mod_v4 as fa
import networkx as nx
from networkx.algorithms.flow import shortest_augmenting_path
from networkx.algorithms.flow import edmonds_karp
import numpy as np
import pickle


bond_file = './set6_x6_eq_bonds.p'
coord_file = './set6_x6_eq_coords.p'

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
for delta in np.arange(0.0,31.0,1.0):
    print 'delta = %2.3f' % delta
    sim.update_bond_coord_lists(direction, delta)
    sim.generate_directed_graph()
    print 'running cut algorithm...'
    cut_value, edge_cut_list = sim.min_cut(direction, attribute, delta, algorithm)
    print 'min cut completed...'
    cut_values.append([delta, cut_value])
    height = sim.find_fracture_height(edge_cut_list, direction)
    print height
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

print "run complete"
pickle.dump(cut_values, open('./cross_link_test/cutvalue_ap-x.p','wb'))
pickle.dump(cut_height, open('./cross_link_test/fractureheight_ap-x.p','wb'))
pickle.dump(cut_freq, open('./cross_link_test/edgefreq_ap-x.p','wb'))
pickle.dump(cut_edges_all, open('./cross_link_test/cutedges_ap-x.p','wb'))