import os
import pandas as pd
import warnings
from pymatgen.io.cif import CifParser
import re

warnings.simplefilter(action="ignore", category=FutureWarning)

def filter_cif_properties(cif_file_path, properties_to_check):
    file_name = os.path.splitext(os.path.basename(cif_file_path))[0]

    parser = CifParser(cif_file_path)
    cif_data = parser.as_dict()

    df = pd.DataFrame(cif_data)

    if not all(property in df[file_name].index for property in properties_to_check):
        return None

    filtered_df = df[file_name][properties_to_check]

    for column_name, column_data in filtered_df.iteritems():
        if column_name in ['_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z']:
            try:
                column_data = [float(re.sub(r'\(\d+\)', '', value)) for value in column_data]
            except ValueError as e:
                print("Error converting value to float:", column_name, column_data)
            filtered_df[column_name] = column_data

    return pd.DataFrame(filtered_df).T

properties_to_check = [
    "_chemical_formula_moiety",
    "_atom_site_fract_x",
    "_atom_site_fract_y",
    "_atom_site_fract_z",
    "_atom_site_type_symbol",
    "_cell_measurement_temperature",
    "_cell_measurement_theta_max",
    "_cell_measurement_theta_min",
    "_cell_volume",
    "_diffrn_radiation_wavelength",
    "_diffrn_reflns_av_R_equivalents",
    "_diffrn_reflns_limit_h_max",
    "_diffrn_reflns_limit_h_min",
    "_diffrn_reflns_limit_k_max",
    "_diffrn_reflns_limit_k_min",
    "_diffrn_reflns_limit_l_max",
    "_diffrn_reflns_limit_l_min",
    "_diffrn_reflns_number",
    "_diffrn_reflns_theta_max",
    "_diffrn_reflns_theta_min",
    "_exptl_absorpt_coefficient_mu",
    "_exptl_absorpt_correction_T_max",
    "_exptl_absorpt_correction_T_min",
    "_exptl_crystal_size_max",
    "_exptl_crystal_size_mid",
    "_exptl_crystal_size_min",
    "_refine_diff_density_max",
    "_refine_diff_density_min",
    "_refine_ls_extinction_coef"
]

folder_path = "/content/Cod_final"
output_csv_file = "/content/CIF_DATA.csv"

cif_files = os.listdir(folder_path)
cif_files = [os.path.join(folder_path, file) for file in cif_files if file.endswith(".cif")]

output_dataframe = pd.DataFrame(columns=["File Name"] + properties_to_check)

for cif_file_path in cif_files:
    filtered_dataframe = filter_cif_properties(cif_file_path, properties_to_check)

    if filtered_dataframe is not None:
        filtered_dataframe["File Name"] = os.path.basename(cif_file_path)

        filtered_dataframe = filtered_dataframe[["File Name"] + properties_to_check]

        output_dataframe = output_dataframe.append(filtered_dataframe, ignore_index=True)
output_dataframe.to_csv(output_csv_file, index=False)

print("CSV file created successfully.")
