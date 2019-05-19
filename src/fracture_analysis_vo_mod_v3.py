from pyparsing import Optional,Combine,Suppress,Group,OneOrMore,Literal,Word,nums
import networkx as nx
import numpy as np
import pickle

class FractureAnalysis():

    def __init__(self, bond_file, coord_file):
        self.g = nx.DiGraph()
        self.reset_graph()
        self.create_bond_dict()
        self.get_coord_info(coord_file)
        self.get_bond_info(bond_file)
        self.updated_bond_list = []
        self.updated_coord_list = []


    def reset_graph(self):
        self.g.clear()

    def create_bond_dict(self):
        self.bond_list = {}
        self.bond_list['1-1'] = 3.8 # C-C
        self.bond_list['1-2'] = 3.8 # C-C
        self.bond_list['1-10'] = 3.8 # C-C
        self.bond_list['1-11'] = 3.8 # C-C
        self.bond_list['2-2'] = 3.8 # C-C
        self.bond_list['2-10'] = 3.8 # C-C
        self.bond_list['10-10'] = 3.8 # C-C
        self.bond_list['2-11'] = 3.8 # C-C
        self.bond_list['10-11'] = 3.8 # C-C
        self.bond_list['11-11'] = 3.8 # C-C
        self.bond_list['1-8'] = 3.7 # C-S
        self.bond_list['2-8'] = 3.7 # C-S
        self.bond_list['8-10'] = 3.7 # C-S
        self.bond_list['8-11'] = 3.7 # C-S
        self.bond_list['1-4'] = 4.0 # C-O
        self.bond_list['1-5'] = 4.0 # C-O
        self.bond_list['1-6'] = 4.0 # C-O
        self.bond_list['2-4'] = 4.0 # C-O
        self.bond_list['2-5'] = 4.0 # C-O
        self.bond_list['2-6'] = 4.0 # C-O
        self.bond_list['4-10'] = 4.0 # C-O
        self.bond_list['4-11'] = 4.0 # C-O
        self.bond_list['5-10'] = 4.0 # C-O
        self.bond_list['5-11'] = 4.0 # C-O
        self.bond_list['6-10'] = 4.0 # C-O
        self.bond_list['6-11'] = 4.0 # C-O
        self.bond_list['1-3'] = 4.0 # C-N
        self.bond_list['1-9'] = 4.0 # C-N
        self.bond_list['1-12'] = 4.0 # C-N
        self.bond_list['2-3'] = 4.0 # C-N
        self.bond_list['2-9'] = 4.0 # C-N
        self.bond_list['2-12'] = 4.0 # C-N
        self.bond_list['3-10'] = 4.0 # C-N
        self.bond_list['3-11'] = 4.0 # C-N
        self.bond_list['9-10'] = 4.0 # C-N
        self.bond_list['10-12'] = 4.0 # C-N
        self.bond_list['11-12'] = 4.0 # C-N
        self.bond_list['9-11'] = 4.0 # C-N
        self.bond_list['4-8'] = 5.0 # S-O
        self.bond_list['5-8'] = 5.0 # S-O
        self.bond_list['6-8'] = 5.0 # S-O

#    def get_bond_info(self, bond_file):
#        # Extracts edge data, including various stress energy data
#        # Note that average stress energy based on the arithmetic average is arbitrarily selected
#        # ---INPUTS---
#        # bond_file: edge data file location
#        
#        #Initialize list that will hold all edge data    
#        bond_info_list = []
#        
#        # Establish grammar
#        bond_id = Word(nums)
#        bond_type = Word(nums)
##        bonded_atom_id= Word(nums)
##        bonded_order= Word(nums)
##        arithmetic_avg_stress = Word(nums)
##        absolute_geometric_avg_stress = Word(nums)
#        bond_distance = Word(nums)
#        lparan = Suppress(Literal('('))
#        rparan = Suppress(Literal(')'))
#        comma = Suppress(Literal(','))
##        num = Word(nums)
#        #floating = Combine(Optional(Literal('0'))+Literal('.')+Word(nums))
#        floating = Combine(Optional(Word(nums))+Literal('.')+Word(nums))
##        bond_info = Group(lparan + Word(nums) + floating + floating + floating + floating +rparan+Optional(comma))
#        bond_info = Group(lparan + Word(nums) + floating + floating + rparan+Optional(comma))
#        bond_parser = bond_id+bond_type+OneOrMore(bond_info)
#        
#        # Parse file
#        bond_fp=open(bond_file)
#        bond_fp.readline()
#        for linestr in bond_fp:
#            try:
#                ## run each line through the parser...
#                result = bond_parser.parseString(linestr)
#                c = [i for i in result]
#                bond_id_val = c[0]
#                bond_id_type = c[1]
#                for val in c[2:]:
#                    attr_dict = {}
#    #                    connecting_atom, bo, arith_avg_stress, geo_avg_stress, bond_distance = tuple(val)
#                    connecting_atom, bo, bond_distance = tuple(val)
#                    #attr_dict['bond_type'] = bond_id_type
#                    attr_dict['bond_distance'] = float(bond_distance)
#                    #attr_dict['arith_avg_stress'] = float(arith_avg_stress)
#                    bond = [int(bond_id_type), int(self.full_coord_info_list[int(connecting_atom)-1][1]['bond_type'])]                   
#                    if bond[0] < bond[1]:
#                        bond_idx = str(bond[0]) + '-' + str(bond[1])
#                    else:
#                        bond_idx = str(bond[1]) + '-' + str(bond[0])
#                    attr_dict['rupture_force'] = self.bond_list[bond_idx]
##                    if 1 in bond or 2 in bond or 10 in bond:
##                        if 3 in bond:
##                            attr_dict['rupture_force'] = 2.8
##                        elif 2 in bond:
##                            attr_dict['rupture_force'] = 3.3
##                    elif bond[0] == 3 and bond[1] == 3:
##                        attr_dict['rupture_force'] = 4.1
##                    else:
##                        attr_dict['rupture_force'] = 0.0
##                        print 'Unsupported bond encountered...'
#                    #attr_dict['rupture_force'] = 3.3
#                    bond_info_list.append([bond_id_val,connecting_atom,attr_dict])
#            except:
#                print 'ERROR'
##                print bond_id_val, connecting_atom, bond
##                break
#        self.full_bond_info_list = bond_info_list

    def get_bond_info(self, bond_file):
        bond_info_list = []

        fp = pickle.load(open(bond_file, 'rb'))
        for linestr in fp:
            try:
                attr_dict = {}
                bond_id1 = linestr[2]
                bond_id2 = linestr[3]
                atom_type1 = int(self.full_coord_info_list[bond_id1-1][1]['bond_type'])
                atom_type2 = int(self.full_coord_info_list[bond_id2-1][1]['bond_type'])
                if (atom_type1 != 7) and (atom_type2 != 7):
                    if atom_type1 < atom_type2:
                        attr_dict['rupture_force'] = self.bond_list[str(atom_type1) + '-' + str(atom_type2)]
                    else:
                        attr_dict['rupture_force'] = self.bond_list[str(atom_type2) + '-' + str(atom_type1)]
                    bond_info_list.append( [bond_id1, bond_id2, attr_dict] )
            except:
                print 'ERROR'
                print bond_id1, bond_id2
                print atom_type1, atom_type2
                attr_dict['rupture_force'] = 10.0
                #break
        self.full_bond_info_list = bond_info_list


    def get_coord_info(self, coord_file):
        # Extract cartesian coordinate data for each node
        # ---INPUTS---
        # coord_file: coordinate file location

        # Initialize list that will hold node data
        attr_list = []

        # Establish grammar
        #ID = Word(nums)
        #mol_ID = Word(nums)
        #Type = Word(nums)
#        x = Word(nums)
#        y = Word(nums)
#        z = Word(nums)
        #integer = Combine(Optional(Literal('-')) + Word(nums))
        #floating = Combine(integer+Literal('.')+Word(nums))

        #bond_info = Group(floating + floating + floating + floating + integer + integer + integer)
        #coord_parser = ID + mol_ID + Type + floating + floating + floating + floating + integer + integer + integer#OneOrMore(bond_info)

        fp=pickle.load(open(coord_file,'rb'))
        count = 1
        for linestr in fp:
            try:
                attr_dict = {}
                #result = coord_parser.parseString(linestr)
                #bond_id_val = result[0]
                bond_id_val = linestr[0]
                if bond_id_val == count:                    
                    #bond_id_type = result[2]
                    bond_id_type = linestr[1]
                    #charge, x_val, y_val, z_val, bx, by, bz = result[3]
                    x_val = linestr[2]
                    y_val = linestr[3]
                    z_val = linestr[4]
                    attr_dict['bond_type'] = bond_id_type
                    attr_dict['x_coord'] = float(x_val)
                    attr_dict['y_coord'] = float(y_val)
                    attr_dict['z_coord'] = float(z_val)
                    attr_list.append([bond_id_val, attr_dict])
                    count = count + 1
                elif bond_id_val > count:
                    while bond_id_val > count:
                        attr_list.append([count, 0])
                        count = count + 1
                    bond_id_type = linestr[1]
                    #charge, x_val, y_val, z_val, bx, by, bz = result[3]
                    x_val = linestr[2]
                    y_val = linestr[3]
                    z_val = linestr[4]
                    attr_dict['bond_type'] = bond_id_type
                    attr_dict['x_coord'] = float(x_val)
                    attr_dict['y_coord'] = float(y_val)
                    attr_dict['z_coord'] = float(z_val)
                    attr_list.append([bond_id_val, attr_dict])
                    count = count + 1                    
            except:
                print 'Error'
                break
        attr_list = self.adjust_coord_bounds(attr_list)
        self.full_coord_info_list = attr_list


    def adjust_coord_bounds(self, coord_list):
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
            if elem[1] == 0:
                continue
            else:
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
            if elem[1] != 0:
                elem[1]['x_coord'] = elem[1]['x_coord'] - x_min
                elem[1]['y_coord'] = elem[1]['y_coord'] - y_min
                elem[1]['z_coord'] = elem[1]['z_coord'] - z_min


        return coord_list


    def update_bond_coord_lists(self, direction, delta):
        self.updated_bond_list = []
        self.updated_coord_list = []
        # Update coordinate list
        included_nodes = []
        for elem in self.full_coord_info_list:
            if elem[1] != 0:
                node = elem[0]
                attr_dict = elem[1]
                if attr_dict[direction] <= self.coord_center + delta and attr_dict[direction] >= self.coord_center - delta:
                    self.updated_coord_list.append([node, attr_dict])
                    included_nodes.append(node)
        # Update bond info list
        for elem in self.full_bond_info_list:
            if elem[1] != 0:
                node1 = elem[0]
                node2 = elem[1]
                attr_dict = elem[2]
                if node1 in included_nodes and node2 in included_nodes:
                    self.updated_bond_list.append([node1, node2, attr_dict])


    def generate_directed_graph(self, remove_uncoordinated=False):
        # Generates a directed graph based on two input files
        # ---INPUTS---
        # bond_file: edge data associated with bond stress energy
        # coord_file: node data associating atom with cartesian coordinates

        if remove_uncoordinated:
            print 'Removing uncoordinated O...'
            self.remove_uncoordinated_atoms('2')

        for elem in self.updated_bond_list:
            bond_id = elem[0]
            connecting_atom = elem[1]
            attr_dict = elem[2]
            self.g.add_edge(bond_id, connecting_atom, attr_dict)

        for elem in self.updated_coord_list:
            node = elem[0]
            attr_dict = elem[1]
            self.g.node[node] = attr_dict



    def find_center(self, direction):
        # Finds the center of the material based on coordinate
        # This function is used in the 'run_max_flow' function to aid in the clamping of nodes
        # outside the specified 'delta' distance
        # ---INPUTS---
        # direction: the coordinate direction of interest, e.g., if the fracture analysis
        # takes place along the x-direction, then direction = 'x_coord'
        coord_dict = []
        for elem in self.full_coord_info_list:
            if elem[1] != 0:
                attr_dict = elem[1]
                coord_dict.append(attr_dict[direction])
                self.max_val = np.max(coord_dict)
                self.min_val = np.min(coord_dict)
        self.coord_center = 0.5*(self.max_val + self.min_val)


    def adjust_graph_capacities(self, direction, attribute, max_edge_value):
        # Adjusts the capacities for edges oriented from sink to source so that
        # they are equal to 2*2*max_edge_value
        for (node1,node2) in self.g.edges():
            if ((node1 == 8307) or (node2 == 8307)):
                print "break"
            node1_pos = nx.get_node_attributes(self.g,direction)[node1]
            node2_pos = nx.get_node_attributes(self.g,direction)[node2]
            periodic_box_length = (self.max_val - self.min_val)/2
            pos_diff = abs(node1_pos - node2_pos)
            hi_pos1_diff = abs(node1_pos - self.max_val)
            hi_pos2_diff = abs(node2_pos - self.max_val)
            lo_pos1_diff = abs(node1_pos - self.min_val)
            lo_pos2_diff = abs(node2_pos - self.min_val)
            if pos_diff < periodic_box_length:
                if node1_pos < node2_pos:
                    self.g.edge[node1][node2][attribute] = 2*2*max_edge_value
            #elif ((hi_pos1_diff < 7.0) or (hi_pos2_diff < 7.0)):
            #    self.g.edge[node1][node2][attribute] = 2*2*max_edge_value
            #elif ((lo_pos1_diff < 7.0) or (lo_pos2_diff < 7.0)):
            #    self.g.edge[node1][node2][attribute] = 2*2*max_edge_value            
            else: #this edge probably crosses a periodic boundary
                self.g.edge[node1][node2][attribute] = 2*2*max_edge_value


    def remove_uncoordinated_atoms(self, atom_idx, debug=True):
        if debug:
            print 'number of nodes before pruning: %d' % len(self.g.nodes())

#        for node in self.g.nodes():
#            bond_type = nx.get_node_attributes(self.g,'bond_type')[node]
#            if bond_type == atom_idx:
#                try:
#                    num_edges = len(self.g.edge[node])
#                except:
#                    num_edges = 0
#                if num_edges < 2:
#                    self.g.remove_node(node)
        for (idx,elem) in enumerate(self.updated_coord_list):
            node = elem[0]
            attr_dict = elem[1]
            bond_type = attr_dict['bond_type']
            if bond_type == atom_idx:
                num_edges = self.count_edges(node, self.updated_bond_list)
                if num_edges < 2:
                    self.updated_coord_list.__delitem__(idx)
                    #print 'removing node: %s' % node
                    if num_edges > 0:
                        for (idx2, elem2) in enumerate(self.updated_bond_list):
                            node1 = elem2[0]
                            node2 = elem2[1]
                            #print node1
                            if node1 == node or node2 == node:
                                #print 'removing associated bond'
                                self.updated_bond_list.__delitem__(idx2)
        if debug:
            print 'number of nodes after pruning: %d' % len(self.g.nodes())


    def count_edges(self, node, bond_list):
        count = 0
        for elem in bond_list:
            node1 = elem[0]
            node2 = elem[1]
            if node1 == node or node2 == node:
                count += 1
        return count


    def add_source_sink_nodes(self):
        self.source = str(len(self.full_coord_info_list) + 1)
        self.sink = str(len(self.full_coord_info_list) + 2)
        self.g.add_node(self.source)
        self.g.add_node(self.sink)


#    def connect_source_sink_to_surfaces(self, direction, attribute, max_stress):
#        top_surface_nodes = []
#        bottom_surface_nodes = []
#        for target in self.g.nodes():
#            if target != self.source or target != self.sink:

#                

    def build_source_sink_list(self, direction):
        # find surface nodes
        top_surface_nodes = []
        bottom_surface_nodes = []
        for target in self.g.nodes():
            top_flag = 1
            bottom_flag = 1
            try:
                for neigh in self.g.neighbors(target):
                    if self.g.node[target][direction] < self.g.node[neigh][direction]:
                        top_flag = 0
                    if self.g.node[target][direction] > self.g.node[neigh][direction]:
                        bottom_flag = 0
                if top_flag == 1:
                    top_surface_nodes.append(target)
                if bottom_flag == 1:
                    bottom_surface_nodes.append(target)
            except:
                pass
        return top_surface_nodes, bottom_surface_nodes

    def build_source_sink_list2(self, direction):
        # find surface nodes
        arr_len = len(self.full_coord_info_list)
        top_coord_array = np.zeros(arr_len)
        bottom_coord_array = np.zeros(arr_len)
        top_surface_nodes = np.zeros(arr_len)
        bottom_surface_nodes = np.zeros(arr_len)
        for target in self.g.nodes():
            top_flag = 1
            bottom_flag = 1
            try:
                for neigh in self.g.neighbors(target):
                    if self.g.node[target][direction] < self.g.node[neigh][direction]:
                        top_flag = 0
                    if self.g.node[target][direction] > self.g.node[neigh][direction]:
                        bottom_flag = 0
                if top_flag == 1:
                    #top_surface_nodes.append(target)
                    top_surface_nodes[int(target)] = target
                    top_coord_array[int(target)]= self.g.node[target][direction]
                if bottom_flag == 1:
                    #bottom_surface_nodes.append(target)
                    bottom_surface_nodes[int(target)] = target
                    bottom_coord_array[int(target)] = self.g.node[target][direction]
                    #bottom_coord_list.append(self.g.node[target][direction])
            except:
                pass

        #top_array = np.asarray(top_coord_list)
        thist,tbedges = np.histogram(top_coord_array,bins=20)
        #bot_array = np.asarray(bottom_coord_list)
        bhist,bbedges = np.histogram(bottom_coord_array,bins=20)
        blo_list = np.where(bhist[1:19]>0)
        blo_ind1 = (blo_list[0][0]) + 1
        blo_ind2 = (blo_list[0][2]) + 1
        
        top_lim = float(tbedges[17])
        top2_lim = float(tbedges[19])
        #bot_lim = float(bbedges[2])
        #bot2_lim = float(bbedges[1])
        bot_lim = float(bbedges[blo_ind1])
        bot2_lim = float(bbedges[blo_ind2])
        lim_top_surface_nodes_list = []
        lim_bottom_surface_nodes_list = []
        for ii in xrange(len(top_coord_array)):
            if ((top_coord_array[ii] > top_lim) and (top_coord_array[ii] < top2_lim)):
                lim_top_surface_nodes_list.append(int(top_surface_nodes[ii]))
        for jj in xrange(len(bottom_coord_array)):
            if bottom_coord_array[jj] != 0.0:
                if ((bottom_coord_array[jj] < bot2_lim) and (bottom_coord_array[jj] > bot_lim)):
                    lim_bottom_surface_nodes_list.append(int(bottom_surface_nodes[jj]))    
        return lim_top_surface_nodes_list, lim_bottom_surface_nodes_list


    def connect_source_sink_to_surfaces(self, top_surface_nodes, bottom_surface_nodes, attribute, max_stress):
        # Link surface nodes to source or sink
        for node in top_surface_nodes:
            self.g.add_edge(self.source,node,{attribute: 2*2*max_stress})
        for node in bottom_surface_nodes:
            self.g.add_edge(node,self.sink,{attribute: 2*2*max_stress})

    def write_nodes_coords(self,top,bottom,edge_cut_list,full_coord_info_list):
        atop = np.asarray(top)
        bot = np.asarray(bottom)
        edge_array = np.asarray(edge_cut_list)
        print len(top)
        print len(bottom)
        print len(edge_cut_list)
        print len(full_coord_info_list)
        atype = 0
        for i in xrange(len(full_coord_info_list)):
            #print full_coord_info_list[i]
            atype = 0
            top_ans = np.where(atop==full_coord_info_list[i][0])
            if len(top_ans[0]) > 0:
                #print top[top_ans[0][0]],top_ans[0][0]
                atype = atype + 4
                #atype = 4
                #print "this is top:", count, full_coord_info_list[i]

            bot_ans = np.where(bot==full_coord_info_list[i][0])
            if len(bot_ans[0]) > 0:
                #print bot[bot_ans[0][0]],bot_ans[0][0]
                #print "this is bot:", count, full_coord_info_list[i]
                atype = atype + 5
                #atype = 5

            ea = np.where(edge_array==full_coord_info_list[i][0])
            if len(ea[0]) > 0:
                #print edge_array[ea[0][0]][ea[1][0]],ea[0][0],ea[1][0]
                #print "this is edge:", count, full_coord_info_list[i]
                atype = atype + 3
                #atype = 3
            if atype == 0:
                atype = int(full_coord_info_list[i][1]['bond_type'])
            print int(full_coord_info_list[i][0]),atype,float(full_coord_info_list[i][1]['x_coord']),float(full_coord_info_list[i][1]['y_coord']),float(full_coord_info_list[i][1]['z_coord'])
        print "done"

    def min_cut(self, direction, attribute, delta, algorithm):
        # Runs max flow algorithm utilizing networkx's minimum cut algorithm
        # This function does the follow:
        # 1. Creates an augmented graph where the thickness of the substrate is determined delta
        # 2. Source and sink nodes are then added to the augmented graph and connected to top and bottom surface nodes
        # 3. A minimum cut algorithm is implemented from the source to sink nodes using the average stress
        # ---INPUTS---
        # direction: direction along which the minimum cut will be performed (note: must be a string: 'x_coord', 'y_coord', or 'z_coord') 
        # delta: height cutoff (note: substrate thickness is given by 2*delta)

        # Find max arith avg stress
        try:
            max_stress = np.array(nx.get_edge_attributes(self.g,attribute).values()).max()
        except:
            max_stress = 0

        # Set capacities for nodes oriented from sink to source
        self.adjust_graph_capacities(direction, attribute, max_stress)

        # build source/sink lists
#        top, bottom = self.build_source_sink_list(direction)
        top, bottom = self.build_source_sink_list2(direction) # this function is similar to the original but with mods to reduce the geometrical location of the surface nodes
        # Add source and sink
        self.add_source_sink_nodes()

        # Connect source and sink
        self.connect_source_sink_to_surfaces(top, bottom, attribute, max_stress)

        # Run max flow starting from created source
        cut_value, partitions = nx.minimum_cut(self.g,s=self.source,t=self.sink,capacity=attribute,flow_func=algorithm) #flow_func=shortest_augmenting_path)

        # Determine edges in cut, exclude source/sink connection
        cut = 0
        edge_cut_list = []
        for p1_node in partitions[0]:
            for p2_node in partitions[1]:
                if self.g.has_edge(p1_node,p2_node):
                    if p1_node != self.source and p1_node != self.sink and p2_node != self.source and p2_node != self.sink:
                        edge_cut_list.append((p1_node,p2_node))
                        cut += nx.get_edge_attributes(self.g,attribute)[(p1_node,p2_node)]

        #self.write_nodes_coords(top,bottom,edge_cut_list,self.full_coord_info_list)
        return cut, edge_cut_list


    def find_fracture_height(self, edge_list, dim):
        coord_vals = []
        for (p1,p2) in edge_list:
            coord_vals.append(nx.get_node_attributes(self.g,dim)[p1])
            coord_vals.append(nx.get_node_attributes(self.g,dim)[p2])
        coord_vals = np.array(coord_vals)
        try:
            coord_max = coord_vals.max()
            coord_min = coord_vals.min()
        except:
            coord_max = 0
            coord_min = 0
        return coord_max - coord_min