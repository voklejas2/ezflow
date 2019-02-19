import pickle as pk



coords_file = "coords_sorted.txt"
cf = open(coords_file,'r')
coords = cf.readlines()

coords_data = []
for line in coords:
    line_list = line.split(' ')
    line_data = []
    for elem in line_list:
        line_data.append(round(float(elem),5))
    coords_data.append(line_data)

pk.dump(coords_data, open('coords_crosslink1.p','wb'))


bonds_file = "bonds.txt"
cf = open(bonds_file,'r')
bonds = cf.readlines()

bonds_data = []
for line in bonds:
    line_list = line.split(' ')
    line_data = []
    for elem in line_list:
        line_data.append(int(elem))
    bonds_data.append(line_data)

pk.dump(bonds_data, open('bonds_crosslink1.p','wb'))