import networkx as nx
import matplotlib.pyplot as plt
from mendeleev import element
from networkx.drawing.layout import fruchterman_reingold_layout
import traceback

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
def create_molecular_graph(file_path, row_index):
  # try:
    data = pd.read_csv(file_path)
    label_to_atom = {}
    element_labels = [
        label.replace("\\", "").replace("'", "").replace('"', "").strip()
        for label in data.loc[row_index, '_atom_site_label'].replace('[', '').replace(']', '').split(', ')
    ]

    elements = [
        element.replace("\\", "").replace("'", "").replace('"', "").strip()
        for element in data.loc[row_index, '_atom_site_type_symbol'].replace('[', '').replace(']', '').split(', ')
    ]

    for i, label in enumerate(element_labels):
        label_to_atom[label] = elements[i]

    # Create an empty graph
    G = nx.Graph()

    # Define column names
    atom_label_columns = [
        '_geom_bond_atom_site_label_1',
        '_geom_bond_atom_site_label_2',
    ]
    bond_length_column = '_geom_bond_distance'

    # Fetch atom labels and bond distances
    atom_labels = {}
    for column in atom_label_columns:
        atom_labels[column] = [
            label.replace("'", "").replace('"', "").strip()
            for label in data.loc[row_index, column].strip("[]'").split(', ')
        ]

    bond_values_str = data.loc[row_index, bond_length_column].replace('[', '').replace(']', '').split(',')
    bond_values = [float(value.strip()) for value in bond_values_str if value]

    # Iterate over each pair of atoms
    for atom_label_1, atom_label_2 in zip(atom_labels['_geom_bond_atom_site_label_1'], atom_labels['_geom_bond_atom_site_label_2']):
        atom_label_1 = atom_label_1.strip()
        atom_label_2 = atom_label_2.strip()

        # Skip adding self-loops
        if atom_label_1 == atom_label_2:
            continue

        # Check if a bond distance is available for the pair
        if (atom_label_1, atom_label_2) in zip(element_labels, element_labels):
            index = list(zip(element_labels, element_labels)).index((atom_label_1, atom_label_2))
            bond_value = bond_values[index]
            G.add_edge(atom_label_1, atom_label_2, distance=bond_value)
        else:
            G.add_edge(atom_label_1, atom_label_2)  # Add edge without bond distance

    # Fetch coordinates
    x_coordinates = data.loc[row_index, '_atom_site_fract_x'].replace('[', '').replace(']', '').split(',')
    y_coordinates = data.loc[row_index, '_atom_site_fract_y'].replace('[', '').replace(']', '').split(',')
    z_coordinates = data.loc[row_index, '_atom_site_fract_z'].replace('[', '').replace(']', '').split(',')

    unique_element_labels = set(element_labels)

    for i, label in enumerate(element_labels):
        x = float(x_coordinates[i].strip())
        y = float(y_coordinates[i].strip())
        z = float(z_coordinates[i].strip())

        # Check if label is a duplicate
        if label in unique_element_labels:
            unique_element_labels.remove(label)  # Remove the label from the set of unique labels
        else:
            continue  # Skip adding new nodes for duplicates

        # Check if node with the same label already exists
        if label in G.nodes():
            # Update the position and atomic number of the existing node
            G.nodes[label]['position'] = (x, y, z)
            G.nodes[label]['atomic_number'] = element_dict[label_to_atom[label]]
        else:
            # Add a new node with the label, position, and atomic number
            G.add_node(label, position=(x, y, z), atomic_number=element_dict[label_to_atom[label]])
    singular_properties = {
        'temperature': data.loc[row_index, '_cell_measurement_temperature'],
        'theta_max': data.loc[row_index, '_cell_measurement_theta_max'],
        'theta_min': data.loc[row_index, '_cell_measurement_theta_min'],
        'cell_volume': data.loc[row_index, '_cell_volume'],
        'radiation_wavelength': data.loc[row_index, '_diffrn_radiation_wavelength'],
        'av_R_equivalents': data.loc[row_index, '_diffrn_reflns_av_R_equivalents'],
        'limit_h_max': data.loc[row_index, '_diffrn_reflns_limit_h_max'],
        'limit_h_min': data.loc[row_index, '_diffrn_reflns_limit_h_min'],
        'limit_k_max': data.loc[row_index, '_diffrn_reflns_limit_k_max'],
        'limit_k_min': data.loc[row_index, '_diffrn_reflns_limit_k_min'],
        'limit_l_max': data.loc[row_index, '_diffrn_reflns_limit_l_max'],
        'limit_l_min': data.loc[row_index, '_diffrn_reflns_limit_l_min'],
        'reflns_number': data.loc[row_index, '_diffrn_reflns_number'],
        'reflns_theta_max': data.loc[row_index, '_diffrn_reflns_theta_max'],
        'reflns_theta_min': data.loc[row_index, '_diffrn_reflns_theta_min'],
        'absorpt_coefficient_mu': data.loc[row_index, '_exptl_absorpt_coefficient_mu'],
        'absorpt_correction_T_max': data.loc[row_index, '_exptl_absorpt_correction_T_max'],
        'absorpt_correction_T_min': data.loc[row_index, '_exptl_absorpt_correction_T_min'],
        'crystal_size_max': data.loc[row_index, '_exptl_crystal_size_max'],
        'crystal_size_mid': data.loc[row_index, '_exptl_crystal_size_mid'],
        'crystal_size_min': data.loc[row_index, '_exptl_crystal_size_min'],
        'diff_density_max': data.loc[row_index, '_refine_diff_density_max'],
        'diff_density_min': data.loc[row_index, '_refine_diff_density_min'],
        'ls_extinction_coef': data.loc[row_index, '_refine_ls_extinction_coef']
    }

    G.graph['properties'] = singular_properties

    # Create the pos dictionary using node coordinates
    pos = nx.spring_layout(G, dim=2)
    return G


 