import pandas as pd
import py3Dmol

# Read the CSV file
data = pd.read_csv("/content/CIF_DATA_modified.csv")

# Extract the coordinates and element symbols from the first row
x_coordinates = data.loc[10, '_atom_site_fract_x'].replace('[', '').replace(']', '').split(',')
y_coordinates = data.loc[10, '_atom_site_fract_y'].replace('[', '').replace(']', '').split(',')
z_coordinates = data.loc[10, '_atom_site_fract_z'].replace('[', '').replace(']', '').split(',')
element_symbols = data.loc[10, '_atom_site_type_symbol'].replace("[", "").replace("]", "").replace("'", "").split(', ')

# Scale up the coordinates by 10x
x_coordinates = [float(x) * 10 for x in x_coordinates]
y_coordinates = [float(y) * 10 for y in y_coordinates]
z_coordinates = [float(z) * 10 for z in z_coordinates]

# Calculate the count of element symbols
element_count = len(element_symbols)

# Create a set of unique elements
unique_elements = set(element_symbols)

# Assign colors to each unique element
element_colors = {
    'H': 'white',
    'C': 'gray',
    'N': 'blue',
    'O': 'red',
    'F': 'green',
    'Cl': 'lime',
    'Br': 'darkorange',
}

# Create the XYZ string representation
xyz_string = f"{element_count}\n"
xyz_string += "Molecule\n"

for i in range(element_count):
    symbol = element_symbols[i]
    x = x_coordinates[i]
    y = y_coordinates[i]
    z = z_coordinates[i]
    xyz_string += f"{symbol} {x} {y} {z}\n"

# Define the visualization function
def visualize_molecule(xyz_string, zoom_level=100):
    view = py3Dmol.view(width=1000, height=800)
    view.addModel(xyz_string, 'xyz')

    # Adjust bond size
    bond_radius = 0.1    # Size of the bond cylinders
    view.setStyle({'stick': {'radius': bond_radius}})

    # Add spheres for atoms
    sphere_radius = 0.4  # Adjust the sphere size here
    for i in range(element_count):
        symbol = element_symbols[i]
        x = x_coordinates[i]
        y = y_coordinates[i]
        z = z_coordinates[i]
        color = element_colors.get(symbol, 'gray')
        view.addSphere({'center': {'x': x, 'y': y, 'z': z}, 'radius': sphere_radius, 'color': color})

    # Use element name as labels for atoms
    for i in range(element_count):
        symbol = element_symbols[i]
        x = x_coordinates[i]
        y = y_coordinates[i]
        z = z_coordinates[i]
        view.addLabel(symbol, {'position': {'x': x, 'y': y, 'z': z}, 'fontColor': 'black', 'showBackground': False})

    view.zoomTo(zoom_level)
    return view.show()

# Visualize the molecule
visualize_molecule(xyz_string, zoom_level=200)
