import numpy as np

dt = np.dtype([('id', 'i4'),('type', 'i2'), ('atom1','i4'),('atom2', 'f4')])
atom_dt = np.dtype([('id', 'i4'),('type', 'i2'),('x', 'f4'), ('y', 'f4'), ('z', 'f4')])

box_size = 20.0
direction = 'x'
bonds_file = "bonds_set6_63.txt"
with open (bonds_file,'r') as bf:
    bonds = bf.readlines()

coords_file = "coords_set6_63.txt"
with open (coords_file,'r') as cf:
    coords = cf.readlines()

num_atoms = 21804
num_bonds = 21873

atoms = np.zeros(num_atoms+1,dtype=atom_dt)

dir_min = 1000.0
dir_max = -1000.0
atom_count = 1
for line in coords:
    sline = line.strip()
    line_list = sline.split(' ')
    aid = int(line_list[0])
    atype = int(line_list[1])
    x = float(line_list[2])
    y = float(line_list[3])
    z = float(line_list[4])
    if aid == atom_count:
        atoms[atom_count]['id'] = aid
        atoms[atom_count]['type'] = atype
        atoms[atom_count]['x'] = x
        atoms[atom_count]['y'] = y 
        atoms[atom_count]['z'] = z
        if dir_min > atoms[atom_count][direction]:
            dir_min = atoms[atom_count][direction]
        if dir_max < atoms[atom_count][direction]:
            dir_max = atoms[atom_count][direction]
        atom_count = atom_count + 1
    elif aid > atom_count:
        while aid > atom_count:
            atoms[atom_count]['id'] = 0
            atom_count = atom_count + 1
        atoms[atom_count]['id'] = aid
        atoms[atom_count]['type'] = atype
        atoms[atom_count]['x'] = x
        atoms[atom_count]['y'] = y 
        atoms[atom_count]['z'] = z
        atom_count = atom_count + 1
        
        

dat1 = np.zeros(num_bonds*3,dtype=dt)
count = 1
countpx = num_bonds + 1
#countnx = 2*(num_bonds) + 1
#countpy = 3*(num_bonds) + 1
#countny = 4*(num_bonds) + 1
#countpxpy = 5*(num_bonds) + 1
#countpxny = 6*(num_bonds) + 1
#countnxny = 7*(num_bonds) + 1
#countnxpy = 8*(num_bonds) + 1    
bonds_data = []
for line in bonds:
    line_list = line.split(' ')
    line_data = []
    b_num = int(line_list[0])
    b_type = int(line_list[1])
    if b_type == 13:
        print "this is a giant xlink"
    b_atom1 = int(line_list[2])
    if b_atom1 == 4516:
        print "tp"
    b_atom2 = int(line_list[3])
    if b_atom2 == 4516:
        print "tp"
    #dat1[count]['id'] = count
    #dat1[count]['type'] = b_type
    #dat1[count]['atom1'] = b_atom1
    #dat1[count]['atom2'] = b_atom2
    #line_data.append([b_num,b_type,b_atom1,b_atom2])
    #bonds_data.append(line_data)
    #count = count + 1    
    
    
    dir_diff = abs(atoms[b_atom1][direction]-atoms[b_atom2][direction])
    if dir_diff < box_size:
        line_data = []
        dat1[count]['id'] = count
        dat1[count]['type'] = b_type
        dat1[count]['atom1'] = b_atom1
        dat1[count]['atom2'] = b_atom2
        line_data.append([b_num,b_type,b_atom1,b_atom2])
        bonds_data.append(line_data)
        count = count + 1      
        
        line_data = []
        bpx_atom1 = b_atom1 + 1*(num_atoms) 
        bpx_atom2 = b_atom2 + 1*(num_atoms)
        dat1[countpx]['id'] = countpx
        dat1[countpx]['type'] = b_type
        dat1[countpx]['atom1'] = bpx_atom1
        dat1[countpx]['atom2'] = bpx_atom2    
        line_data.append([countpx,b_type,bpx_atom1,bpx_atom2])
        countpx = countpx + 1
        bonds_data.append(line_data)        
        
    else: # this is probably a bond that crosses a periodic boundary
        if b_type == 13:
            print "this is a giant xlink"        
        try:
            if (atoms[b_atom1][direction] > atoms[b_atom2][direction]): # atom2 lies toward the lo box boundary
                #b_atom1_id = atoms[b_atom1]['id']
                #new_atom_id = atoms[b_atom2]['id'] + num_atoms
                #if ((abs(atoms[b_atom2]['x'] - dir_min) > 3) and (abs(atoms[b_atom2]['y'] - dir_min) > 3) and (abs(atoms[b_atom2]['z'] - dir_min) > 3)): # this atom isn't too close to the sharp end of the structure
                #if abs(atoms[b_atom2][direction] - dir_min) > 4: # this atom isn't too close to the sharp end of the structure
                    #if abs(atoms[b_atom1][direction] - dir_max) > 4: # this atom isn't too close to the sharp end of the structure
                line_data = []
                bpx_atom1 = b_atom1 
                bpx_atom2 = b_atom2 + 1*(num_atoms)
                dat1[countpx]['id'] = countpx
                dat1[countpx]['type'] = b_type
                dat1[countpx]['atom1'] = bpx_atom1
                dat1[countpx]['atom2'] = bpx_atom2    
                line_data.append([countpx,b_type,bpx_atom1,bpx_atom2])
                countpx = countpx + 1
                bonds_data.append(line_data)                                      
                               
            elif (atoms[b_atom2][direction] > atoms[b_atom1][direction]): # atom1 lies toward the lo box boundary
                #b_atom2_id = atoms[b_atom2]['id']
                #new_atom_id = atoms[b_atom1]['id'] + num_atoms
                #if ((abs(atoms[b_atom1]['x'] - dir_min) > 3) or (abs(atoms[b_atom1]['y'] - dir_min) > 3) or (abs(atoms[b_atom1]['z'] - dir_min) > 3)): # this atom isn't too close to the sharp end of the structure
                #if abs(atoms[b_atom1][direction] - dir_min) > 3:
                    #if abs(atoms[b_atom2][direction] - dir_max) > 3: # this atom isn't too close to the sharp end of the structure
                line_data = []
                bpx_atom1 = b_atom1 + 1*(num_atoms)
                bpx_atom2 = b_atom2 
                dat1[countpx]['id'] = countpx
                dat1[countpx]['type'] = b_type
                dat1[countpx]['atom1'] = bpx_atom1
                dat1[countpx]['atom2'] = bpx_atom2    
                line_data.append([countpx,b_type,bpx_atom1,bpx_atom2])
                countpx = countpx + 1
                bonds_data.append(line_data)                   
                
            #line_data = []
            #bpx_atom1 = b_atom1 + 1*(num_atoms) 
            #bpx_atom2 = b_atom2 + 1*(num_atoms)
            #dat1[countpx]['id'] = countpx
            #dat1[countpx]['type'] = b_type
            #dat1[countpx]['atom1'] = bpx_atom1
            #dat1[countpx]['atom2'] = bpx_atom2    
            #line_data.append([countpx,b_type,bpx_atom1,bpx_atom2])
            #countpx = countpx + 1
            #bonds_data.append(line_data)            
            
        except:
            print "error in finding atoms near PB"
        
            # will need to add a check to see if new (added atoms are already bound to existing atoms at a periodic boundary)
            
    
    #line_data = []
    #bpx_atom1 = b_atom1 + 1*(num_atoms) 
    #bpx_atom2 = b_atom2 + 1*(num_atoms)
    #dat1[countpx]['id'] = countpx
    #dat1[countpx]['type'] = b_type
    #dat1[countpx]['atom1'] = bpx_atom1
    #dat1[countpx]['atom2'] = bpx_atom2    
    #line_data.append([countpx,b_type,bpx_atom1,bpx_atom2])
    #countpx = countpx + 1
    #bonds_data.append(line_data)
    
    #line_data = []
    #bnx_atom1 = b_atom1 + 2*(num_atoms) 
    #bnx_atom2 = b_atom2 + 2*(num_atoms)
    #dat1[countnx]['id'] = countnx
    #dat1[countnx]['type'] = b_type
    #dat1[countnx]['atom1'] = bnx_atom1
    #dat1[countnx]['atom2'] = bnx_atom2        
    #line_data.append([countnx,b_type,bnx_atom1,bnx_atom2])      
    #countnx = countnx + 1
    #bonds_data.append(line_data)
    
    #line_data = []
    #bpy_atom1 = b_atom1 + 3*(num_atoms) 
    #bpy_atom2 = b_atom2 + 3*(num_atoms)
    #dat1[countpy]['id'] = countpy
    #dat1[countpy]['type'] = b_type
    #dat1[countpy]['atom1'] = bpy_atom1
    #dat1[countpy]['atom2'] = bpy_atom2           
    #line_data.append([countpy,b_type,bpy_atom1,bpy_atom2])      
    #countpy = countpy + 1    
    #bonds_data.append(line_data)
 
    #line_data = []
    #bny_atom1 = b_atom1 + 4*(num_atoms) 
    #bny_atom2 = b_atom2 + 4*(num_atoms)
    #dat1[countny]['id'] = countny
    #dat1[countny]['type'] = b_type
    #dat1[countny]['atom1'] = bny_atom1
    #dat1[countny]['atom2'] = bny_atom2      
    #line_data.append([countny,b_type,bny_atom1,bny_atom2])      
    #countny = countny + 1     
    #bonds_data.append(line_data)
    
    #line_data = []
    #bpxpy_atom1 = b_atom1 + 5*(num_atoms) 
    #bpxpy_atom2 = b_atom2 + 5*(num_atoms)
    #dat1[countpxpy]['id'] = countpxpy
    #dat1[countpxpy]['type'] = b_type
    #dat1[countpxpy]['atom1'] = bpxpy_atom1
    #dat1[countpxpy]['atom2'] = bpxpy_atom2          
    #line_data.append([countpxpy,b_type,bpxpy_atom1,bpxpy_atom2])      
    #countpxpy = countpxpy + 1   
    #bonds_data.append(line_data)

    #line_data = []
    #bpxny_atom1 = b_atom1 + 6*(num_atoms) 
    #bpxny_atom2 = b_atom2 + 6*(num_atoms)
    #dat1[countpxny]['id'] = countpxny
    #dat1[countpxny]['type'] = b_type
    #dat1[countpxny]['atom1'] = bpxny_atom1
    #dat1[countpxny]['atom2'] = bpxny_atom2           
    #line_data.append([countpxny,b_type,bpxny_atom1,bpxny_atom2])      
    #countpxny = countpxny + 1        
    #bonds_data.append(line_data)
    
    #line_data = []
    #bnxny_atom1 = b_atom1 + 7*(num_atoms) 
    #bnxny_atom2 = b_atom2 + 7*(num_atoms)
    #dat1[countnxny]['id'] = countnxny
    #dat1[countnxny]['type'] = b_type
    #dat1[countnxny]['atom1'] = bnxny_atom1
    #dat1[countnxny]['atom2'] = bnxny_atom2               
    #line_data.append([countnxny,b_type,bnxny_atom1,bnxny_atom2])      
    #countnxny = countnxny + 1
    #bonds_data.append(line_data)
    
    #line_data = []
    #bnxpy_atom1 = b_atom1 + 8*(num_atoms) 
    #bnxpy_atom2 = b_atom2 + 8*(num_atoms)
    #dat1[countnxpy]['id'] = countnxpy
    #dat1[countnxpy]['type'] = b_type
    #dat1[countnxpy]['atom1'] = bnxpy_atom1
    #dat1[countnxpy]['atom2'] = bnxpy_atom2          
    #line_data.append([countnxpy,b_type,bnxpy_atom1,bnxpy_atom2])          
    #countnxpy = countnxpy + 1            
    #bonds_data.append(line_data)
    

with open ('bonds_set6b4.txt','w') as of:
    count = 1
    for i in xrange(len(dat1)):
        if dat1[i][0] != 0:
            print count,dat1[i][1],dat1[i][2],dat1[i][3]
            print dat1[i][0],dat1[i][1],dat1[i][2],dat1[i][3]
            of.write("%d %d %d %d \n" %(count,dat1[i][1],dat1[i][2],dat1[i][3]))
            count = count + 1
    print "not done yet"