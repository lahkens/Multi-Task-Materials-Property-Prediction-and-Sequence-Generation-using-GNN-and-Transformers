import pandas as pd
import py3Dmol

def visualize_molecule(file_path, row_index):
    # Read the CSV data
    data = pd.read_csv(file_path)

    # Initialize dictionaries for bond distances, coordinates, and label to atom mapping
    bond_distances = {}
    coordinates = {}
    label_to_atom = {}

    # Define column names
    atom_label_columns = [
        '_geom_bond_atom_site_label_1',
        '_geom_bond_atom_site_label_2',
    ]
    bond_length_column = '_geom_bond_distance'
    x_column = '_atom_site_fract_x'
    y_column = '_atom_site_fract_y'
    z_column = '_atom_site_fract_z'
    label_column = '_atom_site_label'
    element_column = '_atom_site_type_symbol'

    # Fetch atom labels and bond distances
    atom_labels = {}
    for column in atom_label_columns:
        atom_labels[column] = data.loc[row_index, column].strip("[]'").split(', ')

    bond_values_str = data.loc[row_index, bond_length_column].replace('[', '').replace(']', '').split(',')
    bond_values = [float(value.strip()) for value in bond_values_str if value]
    for j, bond_value in enumerate(bond_values):
        atom_label_1 = atom_labels['_geom_bond_atom_site_label_1'][j].strip()
        atom_label_2 = atom_labels['_geom_bond_atom_site_label_2'][j].strip()
        bond_distances[f'{atom_label_1} {atom_label_2}'] = bond_value

    # Fetch coordinates and label to atom mapping
    x_coordinates = data.loc[row_index, x_column].replace('[', '').replace(']', '').split(',')
    y_coordinates = data.loc[row_index, y_column].replace('[', '').replace(']', '').split(',')
    z_coordinates = data.loc[row_index, z_column].replace('[', '').replace(']', '').split(',')
    element_labels = data.loc[row_index, label_column].replace("[", "").replace("]", "").replace("'", "").split(', ')
    elements = data.loc[row_index, element_column].replace("[", "").replace("]", "").replace("'", "").split(', ')

    for i, label in enumerate(element_labels):
        x = float(x_coordinates[i].strip())
        y = float(y_coordinates[i].strip())
        z = float(z_coordinates[i].strip())
        coordinates[label] = (x, y, z)

        label_to_atom[label] = elements[i]

    # Create a Py3Dmol view object
    view = py3Dmol.view(width=1000, height=800)

    # Add atoms with original coordinates and labels
    sphere_radius = 0.4
    element_colors = {
        'H': 'white',
        'C': 'gray',
        'N': 'blue',
        'O': 'red',
        'F': 'green',
        'Cl': 'lime',
        'Br': 'darkorange',
    }
    for atom_label, coords in coordinates.items():
        atom_label = atom_label.strip("'")
        symbol = label_to_atom[atom_label]
        x, y, z = coords
        x *= 20
        y *= 20
        z *= 20
        color = element_colors.get(symbol, 'gray')
        view.addSphere({'center': {'x': x, 'y': y, 'z': z}, 'radius': sphere_radius, 'color': color})
        view.addLabel(atom_label, {'position': {'x': x, 'y': y, 'z': z}, 'fontSize': 0.2})

    # Add bond distances as cylinders
    bond_radius = 0.1
    for bond_label, distance in bond_distances.items():
        atom_label_1, atom_label_2 = bond_label.split()
        atom_label_1 = atom_label_1.replace("'", "")
        atom_label_2 = atom_label_2.replace("'", "")
        coords_1 = coordinates[atom_label_1]
        coords_2 = coordinates[atom_label_2]
        x1, y1, z1 = coords_1
        x2, y2, z2 = coords_2
        x1 *= 20
        y1 *= 20
        z1 *= 20
        x2 *= 20
        y2 *= 20
        z2 *= 20
        view.addCylinder({'start': {'x': x1, 'y': y1, 'z': z1}, 'end': {'x': x2, 'y': y2, 'z': z2}, 'radius': bond_radius, 'color': 'gray'})

    # Set the style and render the molecule
    view.setStyle({'stick': {}})
    view.zoomTo()
    view.render()

    # Show the molecule
    view.show()