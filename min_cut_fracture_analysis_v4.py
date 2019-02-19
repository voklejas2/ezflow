from pyparsing import Optional,Combine,Suppress,Group,OneOrMore,Literal,alphas,Word,nums
import pdb
import networkx 
import numpy as np
import matplotlib.pyplot as plt
import pickle
from networkx.algorithms.flow import shortest_augmenting_path
from networkx.algorithms.flow import edmonds_karp

def generate_graph(bond_file, coord_file):
    # Generates a directed graph based on two input files
    # ---INPUTS---
    # bond_file: edge data associated with bond stress energy
    # coord_file: node data associating atom with cartesian coordinates
    G = networkx.DiGraph()
    
    bond_info_list = get_bond_info(bond_file)
    for elem in bond_info_list:
        bond_id = elem[0]
        connecting_atom = elem[1]
        attr_dict = elem[2]
        G.add_edge(bond_id, connecting_atom, attr_dict)
    
    coord_info_list = get_coord_info(coord_file)
    for elem in coord_info_list:
        node = elem[0]
        attr_dict = elem[1]
        G.node[node] = attr_dict
    return G

   
def get_bond_info(bond_file):
    # Extracts edge data, including various stress energy data
    # Note that average stress energy based on the arithmetic average is arbitrarily selected
    # ---INPUTS---
    # bond_file: edge data file location
    
    #Initialize list that will hold all edge data    
    bond_info_list = []
    
    # Establish grammar
    bond_id = Word(nums)
    bond_type = Word(nums)
    bonded_atom_id= Word(nums)
    bonded_order= Word(nums)
    arithmetic_avg_stress = Word(nums)
    absolute_geometric_avg_stress = Word(nums)
    bond_distance = Word(nums)
    lparan = Suppress(Literal('('))
    rparan = Suppress(Literal(')'))
    comma = Suppress(Literal(','))
    num = Word(nums)
    #floating = Combine(Optional(Literal('0'))+Literal('.')+Word(nums))
    floating = Combine(Optional(Word(nums))+Literal('.')+Word(nums))
    bond_info = Group(lparan + Word(nums) + floating + floating + floating + floating +rparan+Optional(comma))
    bond_parser = bond_id+bond_type+OneOrMore(bond_info)
    
    # Parse file
    bond_fp=open(bond_file)
    bond_fp.readline()
    for linestr in bond_fp:
        try:
            ## run each line through the parser...
            result = bond_parser.parseString(linestr)
            c = [i for i in result]
            bond_id_val = c[0]
            bond_id_type = c[1]
            for val in c[2:]:
                attr_dict = {}
                connecting_atom, bo, arith_avg_stress, geo_avg_stress, bond_distance = tuple(val)
#                connecting_atom, bo, bond_distance = tuple(val)
                #attr_dict['bond_type'] = bond_id_type
                attr_dict['bond_distance'] = float(bond_distance)
                attr_dict['arith_avg_stress'] = float(arith_avg_stress)
                attr_dict['rupture_force'] = 3.3
                bond_info_list.append([bond_id_val,connecting_atom,attr_dict])
        except:
            print 'ERROR'
            #break
    return bond_info_list

def get_coord_info(coord_file):
    # Extract cartesian coordinate data for each node
    # ---INPUTS---
    # coord_file: coordinate file location

    # Initialize list that will hold node data
    attr_list = []
    
    # Establish grammar
    ID = Word(nums)
    Type = Word(nums)
    x = Word(nums)
    y = Word(nums)
    z = Word(nums)
    integer = Combine(Optional(Literal('-')) + Word(nums))
    floating = Combine(integer+Literal('.')+Word(nums))

    bond_info = Group(floating + floating + floating)
    coord_parser = ID + Type + OneOrMore(bond_info)

    fp=open(coord_file)
    for linestr in fp:
        try:
            attr_dict = {}
            result = coord_parser.parseString(linestr)
            bond_id_val = result[0]
            bond_id_type = result[1]
            x_val, y_val, z_val = result[2]
            attr_dict['bond_type'] = bond_id_type
            attr_dict['x_coord'] = float(x_val)
            attr_dict['y_coord'] = float(y_val)
            attr_dict['z_coord'] = float(z_val)
            attr_list.append([bond_id_val, attr_dict])
        except:
            print 'Error'
            break
    attr_list = adjust_coord_bounds(attr_list)
    return attr_list



def adjust_coord_bounds(coord_list):
    # Adjusts coordinate data so that the minimum x,y,z values are zero
    # This basically done for convenience
    # ---INPUTS---
    # coord_list: list of coordinates, format set by 'get_coord_info' function

    # Initializes coordinate lists
    x_list = []
    y_list = []
    z_list = []
    
    # Locate minimum coordinate values
    for elem in coord_list:
        attr_dict = elem[1]
        x_list.append(attr_dict['x_coord'])
        y_list.append(attr_dict['y_coord'])
        z_list.append(attr_dict['z_coord'])
    x_list = np.array(x_list)
    y_list = np.array(y_list)
    z_list = np.array(z_list)
    x_min = x_list.min()
    y_min = y_list.min()
    z_min = z_list.min()

    # Adjust coordinates so that minimum value is zero
    # elem = [bond_id, attr_dictionary]
    for elem in coord_list:
        elem[1]['x_coord'] = elem[1]['x_coord'] - x_min
        elem[1]['y_coord'] = elem[1]['y_coord'] - y_min
        elem[1]['z_coord'] = elem[1]['z_coord'] - z_min
    return coord_list
    


def find_center(graph,direction):
    # Finds the center of the material based on coordinate
    # This function is used in the 'run_max_flow' function to aid in the clamping of nodes
    # outside the specified 'delta' distance
    # ---INPUTS---
    # graph: the graph to examine
    # direction: the coordinate direction of interest, e.g., if the fracture analysis
    # takes place along the x-direction, then direction = 'x_coord'
    coord_dict = networkx.get_node_attributes(graph,direction)
    coord_list = np.array(coord_dict.values())
    coord_center = 0.5*(coord_list.max()+coord_list.min())
    return coord_center

def adjust_graph_capacities(graph,direction,attribute,max_edge_value):
    # Adjusts the capacities for edges oriented from sink to source so that
    # they are equal to 2*2*max_edge_value
    for (node1,node2) in graph.edges():
        node1_pos = networkx.get_node_attributes(graph,direction)[node1]
        node2_pos = networkx.get_node_attributes(graph,direction)[node2]
        if node1_pos < node2_pos:
            graph.edge[node1][node2][attribute] = 2*2*max_edge_value
    return graph


def find_fracture_height(graph, edge_list, dim):
    coord_vals = []
    for (p1,p2) in edge_list:
        coord_vals.append(networkx.get_node_attributes(graph,dim)[p1])
        coord_vals.append(networkx.get_node_attributes(graph,dim)[p2])
    coord_vals = np.array(coord_vals)
    try:
        coord_max = coord_vals.max()
        coord_min = coord_vals.min()
    except:
        coord_max = 0
        coord_min = 0
    return coord_max - coord_min


def run_max_flow(graph, direction, attribute, delta, algorithm=edmonds_karp):
    # Runs max flow algorithm utilizing networkx's minimum cut algorithm
    # This function does the follow:
    # 1. Creates an augmented graph where the thickness of the substrate is determined delta
    # 2. Source and sink nodes are then added to the augmented graph and connected to top and bottom surface nodes
    # 3. A minimum cut algorithm is implemented from the source to sink nodes using the average stress
    # ---INPUTS---
    # graph: graph describing substrate
    # direction: direction along which the minimum cut will be performed (note: must be a string: 'x_coord', 'y_coord', or 'z_coord') 
    # delta: height cutoff (note: substrate thickness is given by 2*delta)
    
    # Find location of center with respect to selected cartesian coordinate
    center = find_center(graph,direction)
    
    #adjust graph to only include data points within cutoff from center
    aug_graph = graph.subgraph([n for n,attr_dict in graph.node.items() if attr_dict[direction] <= center + delta and attr_dict[direction] >= center - delta])
    vals = np.array(networkx.get_node_attributes(aug_graph,direction).values())
    
    
    # remove oxygens with single edge
    print 'number of nodes before pruning: %d' % len(aug_graph.nodes())    
    for node in aug_graph.nodes():
        bond_type = networkx.get_node_attributes(aug_graph,'bond_type')[node]
        if bond_type == '2':
            num_edges = len(aug_graph.edge[node])
            if num_edges < 2:
                aug_graph.remove_node(node)
    
    print 'number of nodes after pruning: %d' % len(aug_graph.nodes())    
    
    # find surface nodes
    top_surface_nodes = []
    bottom_surface_nodes = []
    for target in aug_graph.nodes():
        top_flag = 1
        bottom_flag = 1
        for neigh in aug_graph.neighbors(target):
            if aug_graph.node[target][direction] < aug_graph.node[neigh][direction]:
                top_flag = 0
            if aug_graph.node[target][direction] > aug_graph.node[neigh][direction]:
                bottom_flag = 0
        if top_flag == 1:
            top_surface_nodes.append(target)
        if bottom_flag == 1:
            bottom_surface_nodes.append(target)
            
    # add source and sink nodes (note that nodes in augmented graph are not relabeled)
    source = str(len(graph.nodes()) + 1)
    sink = str(len(graph.nodes()) + 2)
    aug_graph.add_node(source)
    aug_graph.add_node(sink)
    
    print 'Graph has %d nodes and %d edges...' % (len(aug_graph.nodes()), len(aug_graph.edges()))
#    print 'Source node: %s, Target node: %s' % (source, sink)
    
    # Find max arith avg stress
    try:
        max_stress = np.array(networkx.get_edge_attributes(aug_graph,attribute).values()).max()
    except:
        max_stress = 0
    #print networkx.get_edge_attributes(aug_graph,'arith_avg_stress').values()
    
    # Set capacities for nodes oriented from sink to source
    adjust_graph_capacities(aug_graph,direction,attribute,max_stress)
    
    # Link surface nodes to source or sink
    # edge weight for source/sink is 2*(dim-1)*max_stress value
    for node in top_surface_nodes:
        aug_graph.add_edge(source,node,{attribute: 2*2*max_stress})
    for node in bottom_surface_nodes:
        aug_graph.add_edge(node,sink,{attribute: 2*2*max_stress})

    
    # Run max flow starting from created source
    cut_value, partitions = networkx.minimum_cut(aug_graph,s=source,t=sink,capacity=attribute,flow_func=algorithm) #flow_func=shortest_augmenting_path)
    
    # Determine edges in cut, exclude source/sink connection
    cut = 0
    edge_cut_list = []
    for p1_node in partitions[0]:
        for p2_node in partitions[1]:
            if aug_graph.has_edge(p1_node,p2_node):
                if p1_node != source and p1_node != sink and p2_node != source and p2_node != sink:
                    edge_cut_list.append((p1_node,p2_node))
                    cut += networkx.get_edge_attributes(aug_graph,attribute)[(p1_node,p2_node)]
    
    return cut, edge_cut_list, aug_graph


        


#--------------
# Specifications
#--------------
bond_file = '.\data160128\data160128\center_void_quartz_bond_topology_stress_all.txt'
coord_file = '.\coord_center.txt'
#bond_file = "./data170501/topology_500.txt"
#coord_file = './data170501/500_coords'
flow_dir = 'z_coord'
edge_weight_attr = 'rupture_force' # Value of 3.3 nN is hard-coded into function get_bond_info
algorithm = shortest_augmenting_path  # Can also be edmonds_karp or preflow_push

#---------------
# Generate graph based on molecular structure
#---------------
g = generate_graph(bond_file,coord_file)

#num=0
#for i in range(len(g.edge)):
#    edges = len(g.edge[str(i+1)])
#    if edges < 3:
#        num += 1
#print num
#
#num1 = 0
#num2 = 0
#for node in g.nodes():
#    bond_type = networkx.get_node_attributes(g,'bond_type')[node]
#    if bond_type == '1':
#        num1 += 1
#    else:
#        num2 += 1
#print num1, num2


#---------------
# Run algorithm over a range of Delta values
#---------------
cut_values = []
cut_edges = []
cut_freq = []
cut_height = []
cut_edges_all = []
for delta in np.arange(0.1,3.5,0.02):
    print 'delta = %2.3f' % delta
    cut_value, edge_cut_list, aug_graph = run_max_flow(graph=g,direction=flow_dir,attribute=edge_weight_attr,delta=delta,algorithm=algorithm)
    cut_values.append([delta, cut_value])
    height = find_fracture_height(aug_graph, edge_cut_list, flow_dir)
    cut_height.append([delta, height])
    for edge in edge_cut_list:
        if edge in cut_edges:
            ind = cut_edges.index(edge)
            cut_freq[ind] += 1
        else:
            cut_edges.append(edge)
            cut_freq.append(1)
    cut_edges_all.append([delta,edge_cut_list])
    print ''


#---------------
# Save data
#---------------
#pickle.dump(cut_values, open('center_void_quartz_cutvalue_ek-y.p','wb'))
#pickle.dump(cut_height, open('center_void_quartz_fractureheight_ek-y.p','wb'))
#pickle.dump(cut_freq, open('center_void_quartz_edgefreq_ek-y.p','wb'))
#pickle.dump(cut_edges_all, open('center_void_quartz_cutedges_ek-y.p','wb'))

pickle.dump(cut_values, open('center_void_quartz_cutvalue_ap-z.p','wb'))
pickle.dump(cut_height, open('center_void_quartz_fractureheight_ap-z.p','wb'))
pickle.dump(cut_freq, open('center_void_quartz_edgefreq_ap-z.p','wb'))
pickle.dump(cut_edges_all, open('center_void_quartz_cutedges_ap-z.p','wb'))

