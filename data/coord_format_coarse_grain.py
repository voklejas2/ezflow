import pickle as pk
import numpy as np
import os,sys

dt = np.dtype([('id', 'i4'),('type', 'i2'), ('atom1','i4'),('atom2', 'i4')])
xt = np.dtype([('id', 'i4'),('type', 'i2'), ('atom1','i4'),('atom2', 'i4'),('atom3', 'i4'),('atom4', 'i4')])
atom_dt = np.dtype([('id', 'i4'),('type', 'i2'),('x', 'f4'), ('y', 'f4'), ('z', 'f4')])

def find_neighbors(atom0,bonds):
    neighbor_list = []
    t1 = np.where(bonds['atom1']==atom0)
    t2 = np.where(bonds['atom2']==atom0)
    for nn in xrange(len(t1[0])):
        neighbor_list.append(int(bonds[t1[0][nn]]['atom2']))
    for nn in xrange(len(t2[0])):
        neighbor_list.append(int(bonds[t2[0][nn]]['atom1']))
    
    return(neighbor_list)

def find_xlinks(atom0,atoms,bonds):
    num_links = 0
    link_list = []
    neighbor_atoms = find_neighbors(atom0,bonds)
    # don't count H atoms bonded to N atoms
    legit_atoms = []
    for ila in xrange(len(neighbor_atoms)):
        if atoms[neighbor_atoms[ila]]['type'] != 7:
            legit_atoms.append(neighbor_atoms[ila])
    active_atom = atom0
    chain = []
    chain.append(active_atom)
    while num_links < len(legit_atoms):
        current_neighbors = find_neighbors(active_atom,bonds)
        # add a statement that checks to see if only current neighbors are already in the chain if so then break
        unique_atoms = set(current_neighbors).difference(chain)
        new_atoms = []
        for ua in unique_atoms:
            if ((atoms[ua]['type'] != 4) and (atoms[ua]['type'] != 7)):
                new_atoms.append(ua)
        if len(new_atoms) > 0:            
            for na in current_neighbors:
                if na not in chain:
                #if na == atom0:
                #    continue
                    typen = atoms[na]['type']
                    if typen == 9: # this a connected N atom 
                        link_list.append(na)   
                        chain.append(na)
                        num_links = num_links + 1
                        next
                    elif typen == 12: # this is a connected N atom
                        link_list.append(na)
                        num_links = num_links + 1
                        chain.append(na)
                        next                        
                    elif typen == 11: # this is a connected xlink C atom
                        link_list.append(na)
                        num_links = num_links + 1
                        chain.append(na)
                        next
                    elif ((typen > 1) and (typen != 7) and (typen != 4)):
                        atom0 = active_atom
                        active_atom = na
                        chain.append(na)
                        break
                    elif typen == 1:
                        atom0 = active_atom
                        active_atom = na
                        chain.append(na)
                        continue                   
        else: 
            print "a dead end was reached at xlink atom: ",active_atom
            link_list.append(na)
            num_links = len(legit_atoms)
        
    return(link_list)
    

###########################################
def find_epoxy(atom0,atoms,bonds):
    num_links = 0
    link_list = []
    neighbor_atoms = find_neighbors(atom0,bonds)
    # don't count H atoms bonded to N atoms
    legit_atoms = []
    for ila in xrange(len(neighbor_atoms)):
        if atoms[neighbor_atoms[ila]]['type'] != 7:
            legit_atoms.append(neighbor_atoms[ila])
    active_atom = atom0
    chain = []
    chain.append(active_atom)
    while num_links < len(legit_atoms):
        current_neighbors = find_neighbors(active_atom,bonds)
        # add a statement that checks to see if only current neighbors are already in the chain if so then break
        unique_atoms = set(current_neighbors).difference(chain)
        new_atoms = []
        for ua in unique_atoms:
            if ((atoms[ua]['type'] != 4) and (atoms[ua]['type'] != 7)):
                new_atoms.append(ua)
        if len(new_atoms) > 0:            
            #if active_atom == 13041:
                #print "oh boy"
            for na in current_neighbors:
                if na not in chain:
                #if na == atom0:
                #    continue
                    typen = atoms[na]['type']
                    if typen == 9: # this a connected N atom 
                        #link_list.append(na)   
                        chain.append(na)
                        num_links = num_links + 1
                        next                     
                    elif typen == 12: # this is a connected xlink N atom
                        #link_list.append(na)
                        num_links = num_links + 1
                        chain.append(na)
                        next
                    elif typen == 11: # this is a connected xlink C atom
                        link_list.append(na)
                        num_links = num_links + 1
                        chain.append(na)
                        next       
                    elif typen == 10: # this is a connected xlink C atom
                        link_list.append(na)
                        num_links = num_links + 1
                        chain.append(na)
                        next                         
                    #elif typen == 1:
                    #    atom0 = active_atom
                    #    active_atom = na
                    #    chain.append(na)
                    #    continue                                    
                    elif ((typen > 1) and (typen != 7) and (typen != 4)):
                        nna = find_neighbors(na,bonds)
                        legit_nna = []
                        if typen == 2:
                            for ila in xrange(len(nna)):
                                if atoms[nna[ila]]['id'] != active_atom:
                                    if atoms[nna[ila]]['type'] == 2:
                                        # look into the future
                                        nnna = find_neighbors(nna[ila],bonds)
                                        legit_nnna = []
                                        for nila in xrange(len(nnna)):
                                            if ((atoms[nnna[nila]]['type'] != 7) and (atoms[nnna[nila]]['id'] != na) and (atoms[nnna[nila]]['id'] != atom0)):
                                                legit_nnna.append(nnna[nila])
                                        if len(legit_nnna) > 0:
                                            legit_nna.append(nna[ila])
                                            #atom0 = active_atom
                                            #active_atom = na
                                            #chain.append(na)
                                            #break                 
                                    elif ((atoms[nna[ila]]['type'] != 7) and (atoms[nna[ila]]['type'] != 4) and (atoms[nna[ila]]['id'] != active_atom)):
                                        legit_nna.append(nna[ila])
                        elif typen != 7:
                            for iila in xrange(len(nna)):
                                if ((atoms[nna[iila]]['type'] != 7) and (atoms[nna[iila]]['type'] != 4) and (atoms[nna[iila]]['id'] != active_atom)):                           
                                    legit_nna.append(nna[iila])
                        if len(legit_nna) > 0:
                            atom0 = active_atom
                            active_atom = na
                            chain.append(na)
                            break
                        else:
                            #atom0 = active_atom
                            #active_atom = na
                            chain.append(na)    
                    elif typen == 1:
                        atom0 = active_atom
                        active_atom = na
                        chain.append(na)
                        continue                                       
                    
        else: 
            print "a dead end was reached at epoxy atom: ",active_atom
            link_list.append(na)
            num_links = len(legit_atoms)
        
    return(link_list)
    

###########################################
def main(argv=None):
    num_atoms = 43608
    num_bonds = 43515
    atoms = np.zeros(num_atoms+1,dtype=atom_dt)
    bonds = np.zeros(num_bonds+1,dtype=dt)
    coords_file = "coords_set6_63b2.txt"
    cf = open(coords_file,'r')
    coords = cf.readlines()
    
    count  = 1
    for line in coords:
        sline = line.strip()
        line_list = sline.split(' ')
        coord_id = int(line_list[0])
        coord_type = int(line_list[1])
        x = float(line_list[2])
        y = float(line_list[3])
        z = float(line_list[4])    
        if count == coord_id:
            atoms[count]['id'] = coord_id
            atoms[count]['type'] = coord_type
            atoms[count]['x'] = x
            atoms[count]['y'] = y 
            atoms[count]['z'] = z        
            count = count + 1
        elif count < coord_id:
            while count < coord_id:
                atoms[count]['id'] = count
                atoms[count]['type'] = 0
                atoms[count]['x'] = 0
                atoms[count]['y'] = 0 
                atoms[count]['z'] = 0     
                count = count + 1
            atoms[count]['id'] = coord_id
            atoms[count]['type'] = coord_type
            atoms[count]['x'] = x
            atoms[count]['y'] = y 
            atoms[count]['z'] = z                    
            count = count + 1
        #for elem in line_list:
        #    line_data.append(round(float(elem),5))
        #coords_data.append(line_data)
    
    #pk.dump(coords_data, open('coords_set6_63b2.p','wb'))
    
    
    bonds_file = "bonds_set6b4.txt"
    cf = open(bonds_file,'r')
    bonds_data = cf.readlines()
    
    b_count = 1
    for line in bonds_data:
        sline = line.strip()
        line_list = sline.split(' ')
        bonds[b_count]['id'] = int(line_list[0])
        bonds[b_count]['type'] = int(line_list[1])
        bonds[b_count]['atom1'] = int(line_list[2])
        bonds[b_count]['atom2'] = int(float(line_list[3]))
        b_count = b_count + 1
    
    # now coarse grain it!
    num_xlinks = np.where(bonds['type']==13)
    xlinks = np.zeros(2000,dtype=xt)
    xcount = 1
    for bond in bonds:
        if bond['atom1'] == 2476:
            print "uh-oh"
        if int(bond['type']) == 13:
            #print "now looking at xlink atoms: ",bond['atom1'],bond['atom2']
            # now this is a giant xlink
            # find the neighboring xlink partners
            # find xlink atoms on the other side of the molecule- depending on the type of molecule
            try:
                if atoms[bond['atom1']]['type'] == 12:  # this is a N atom in an xlink molecule
                    neighbor_atoms = find_xlinks(bond['atom1'],atoms,bonds)
                    xlinks[xcount]['id'] = xcount
                    xlinks[xcount]['type'] = int(bond['type'])
                    xlinks[xcount]['atom1'] = int(bond['atom1'])
                    #xlinks[xcount]['atom2'] = int(bond['atom2'])
                    for ix in xrange(len(neighbor_atoms)):
                        new_int = ix + 2
                        tmp_str = "atom" + str(new_int)
                        xlinks[xcount][tmp_str] = neighbor_atoms[ix]                     
                    xcount = xcount + 1
                    #print atoms[bond['atom1']]
                elif atoms[bond['atom1']]['type'] == 9:  # this is a N atom in an xlink molecule
                    neighbor_atoms = find_xlinks(bond['atom1'],atoms,bonds)
                    xlinks[xcount]['id'] = xcount
                    xlinks[xcount]['type'] = int(bond['type'])
                    xlinks[xcount]['atom1'] = int(bond['atom1'])
                    #xlinks[xcount]['atom2'] = int(bond['atom2'])
                    for ix in xrange(len(neighbor_atoms)):
                        new_int = ix + 2
                        tmp_str = "atom" + str(new_int)
                        xlinks[xcount][tmp_str] = neighbor_atoms[ix]                   
                    xcount = xcount + 1                    
                elif atoms[bond['atom2']]['type'] == 12:  # this is a N atom in an xlink molecule
                    neighbor_atoms = find_xlinks(bond['atom2'],atoms,bonds)
                    xlinks[xcount]['id'] = xcount
                    xlinks[xcount]['type'] = int(bond['type'])
                    xlinks[xcount]['atom1'] = int(bond['atom1'])
                    #xlinks[xcount]['atom2'] = int(bond['atom2'])                    
                    for ix in xrange(len(neighbor_atoms)):
                        new_int = ix + 2
                        tmp_str = "atom" + str(new_int)
                        xlinks[xcount][tmp_str] = neighbor_atoms[ix]                   
                    xcount = xcount + 1                    
                elif atoms[bond['atom2']]['type'] == 9:  # this is a N atom in an xlink molecule
                    neighbor_atoms = find_xlinks(bond['atom2'],atoms,bonds)
                    xlinks[xcount]['id'] = xcount
                    xlinks[xcount]['type'] = int(bond['type'])
                    xlinks[xcount]['atom1'] = int(bond['atom1'])
                    #xlinks[xcount]['atom2'] = int(bond['atom2'])                    
                    for ix in xrange(len(neighbor_atoms)):
                        new_int = ix + 2
                        tmp_str = "atom" + str(new_int)
                        xlinks[xcount][tmp_str] = neighbor_atoms[ix]                   
                    xcount = xcount + 1                    
                    
                if atoms[bond['atom1']]['type'] == 11:  # this is a N atom in an xlink molecule
                    neighbor_atoms = find_epoxy(bond['atom2'],atoms,bonds)
                    xlinks[xcount]['id'] = xcount
                    xlinks[xcount]['type'] = int(bond['type'])
                    xlinks[xcount]['atom1'] = int(bond['atom1'])
                    #xlinks[xcount]['atom2'] = int(bond['atom2'])                    
                    for ix in xrange(len(neighbor_atoms)):
                        new_int = ix + 2
                        tmp_str = "atom" + str(new_int)
                        xlinks[xcount][tmp_str] = neighbor_atoms[ix]                   
                    xcount = xcount + 1          
                elif atoms[bond['atom2']]['type'] == 11:  # this is a N atom in an xlink molecule
                    neighbor_atoms = find_epoxy(bond['atom2'],atoms,bonds)
                    xlinks[xcount]['id'] = xcount
                    xlinks[xcount]['type'] = int(bond['type'])
                    xlinks[xcount]['atom1'] = int(bond['atom1'])
                    #xlinks[xcount]['atom2'] = int(bond['atom2'])                    
                    for ix in xrange(len(neighbor_atoms)):
                        new_int = ix + 2
                        tmp_str = "atom" + str(new_int)
                        xlinks[xcount][tmp_str] = neighbor_atoms[ix]                   
                    xcount = xcount + 1                                                  
            except: 
                print "FUCK!",bond
            print "done looking at xlink atoms: ",bond['atom1'],bond['atom2']
    
    count = 1
    for xl in xrange(len(xlinks)):
        for ixl in ['atom2','atom3','atom4']:
            if xlinks[xl][ixl] != 0:
                print xlinks[xl]['atom1']-1,xlinks[xl][ixl]-1
                count = count + 1
    
    print "done?"
    #pk.dump(bonds_data, open('bonds_set6_63b2.p','wb'))

###############################################################################    
if __name__ == "__main__":
    sys.exit( main() )
########################################################################                