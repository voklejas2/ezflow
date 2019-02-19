from __future__ import print_function, division
import numpy as np

atom_dt = np.dtype([('id', 'i4'),('type', 'i2'), ('charge','f4'), ('seg_num','i4'), ('x', 'f4'), ('y', 'f4'), ('z', 'f4')])
# read coord data
f = "coords_glass_all.txt"
node_data = np.zeros(6912,dtype=atom_dt)
with open(f, 'r') as df:
    x = df.readlines()


count = 0
for ia in xrange(6912):
    data_line = x[ia]
    strata_line = data_line.strip()
    if strata_line:  
        spline = data_line.split()
        node_data['type'][count] = int(spline[1])
        node_data['x'][count] = float(spline[2])
        node_data['y'][count] = float(spline[3])
        node_data['z'][count] = float(spline[4])
        node_data['id'][count] = int(spline[0])
        count = count + 1   

        
# read edge data    
print("now read data")
e = "glass_delta_19.6_nodes_lammps.txt"
edge_data = np.zeros(8000,dtype=atom_dt)
with open(e, 'r') as df:
    ex = df.readlines()
    
count = 0
for i in xrange(len(ex)):
    data_line = ex[i]
    ey = ex[i].split(",")
    ez1 = ey[0].strip("()")
    ez2 = ey[1].strip()
    ez3 = ez2.strip(")")
    ez4 = ez3.strip("'")
    eaa = ez1.strip("'")
    ind1 = int(eaa) - 1
    ind2 = int(ez4) - 1
    if node_data[ind1]['type'] == 1:
        node_data[ind1]['type'] = 3
    elif node_data[ind1]['type'] == 2:
        node_data[ind1]['type'] = 4     
    #print(node_data[int(data_line[2:6])]['id'],node_data[int(data_line[2:6])]['type'],node_data[int(data_line[2:6])]['x'],node_data[int(data_line[2:6])]['y'],node_data[int(data_line[2:6])]['z'])
    #print(node_data[int(data_line[10:14])]['id'],node_data[int(data_line[10:14])]['type'],node_data[int(data_line[10:14])]['x'],node_data[int(data_line[10:14])]['y'],node_data[int(data_line[10:14])]['z'])
    count = count + 1

for ii in xrange(len(node_data)):
    print(node_data[ii]['id'],node_data[ii]['type'],node_data[ii]['x'],node_data[ii]['y'],node_data[ii]['z'])