import pandas as pd
import networkx as nx
import numpy as np

element_dict = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
    'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Ni': 27, 'Co': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40,
    'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
    'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
    'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
    'Tl': 81, 'Pb': 82, 'Bi': 83, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96,
    'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106,
    'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116,
    'Ts': 117, 'Og': 118
}

data = pd.read_csv("/content/CIF_DATA_modified.csv")

# Extract the coordinates and element symbols from the first row
x_coordinates = data.loc[4, '_atom_site_fract_x'].replace('[', '').replace(']', '').split(',')
y_coordinates = data.loc[4, '_atom_site_fract_y'].replace('[', '').replace(']', '').split(',')
z_coordinates = data.loc[4, '_atom_site_fract_z'].replace('[', '').replace(']', '').split(',')
element_symbols = data.loc[4, '_atom_site_type_symbol'].replace("[", "").replace("]", "").replace("'", "").split(', ')
print(len(element_symbols))
coordinates = list(zip(x_coordinates, y_coordinates, z_coordinates))

import networkx as nx
import numpy as np

def create_molecular_graph(coordinates, elements):
    num_atoms = len(coordinates)
    atomic_numbers = [element_dict[element] for element in elements]  # Convert element symbols to atomic numbers
    atomic_positions = np.array(coordinates)

    graph = nx.Graph()

    for i in range(num_atoms):
        atomic_number = atomic_numbers[i]
        position = atomic_positions[i]
        graph.add_node(i, atomic_number=atomic_number, position=position)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            graph.add_edge(i, j)

    return graph

molecular_graph = create_molecular_graph(coordinates, element_symbols)

node_attributes = nx.get_node_attributes(molecular_graph, 'atomic_number')
print(node_attributes) 
print(len(node_attributes))

edge_list = list(molecular_graph.edges())
print(edge_list)  
print(len(edge_list))