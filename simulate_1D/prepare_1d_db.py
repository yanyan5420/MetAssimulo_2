import pandas as pd
import numpy as np
import argparse
import json
import pathlib

from read_parameters import read_param
from read_1d_spectra import read_1d_data
from match_names import db_match_cons, input_match_db
from preprocess_1d_spectra import remove_water_calibration, baseline_correction, smooth_spectra, norm_spectra


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--parameter", dest="parameter_file", required=True, help="Input Parameter Files")
args = parser.parse_args()
param_dict = read_param(args.parameter_file)
file_path_1d = param_dict['file_path_1D']
sop_type = param_dict['sop_type']
pulseProgram_type = param_dict['pulseProgrm_type']
# hmdb_dict_path =

base_path = pathlib.Path(__file__).resolve().parents[1]
with open(base_path.joinpath("Input/hmdb_id_names.json")) as json_file:
    hmdb_dict = json.load(json_file)
cons_df_1 = pd.read_csv(base_path.joinpath("Input/cons_df_1.csv"), index_col=0)

data_dict, ppm_scale = read_1d_data(file_path_1d, sop_type, pulseProgram_type)  # read 1d data
print(np.max(data_dict["Citric Acid"]))

match_data_dict = db_match_cons(data_dict, cons_df_1, hmdb_dict)  # format and match db names with cons files
# preprocess 1d data
removed_data_dict = remove_water_calibration(match_data_dict, ppm_scale, [4.5, 5.0], [-0.3, 0.3])
corrected_data_dict = baseline_correction(removed_data_dict, 16, 0.1)
smooth_data_dict = smooth_spectra(corrected_data_dict, 0.05)
norm_data_dict = norm_spectra(smooth_data_dict)

# sample pka for the db
pka_dict = dict()
for meta_name, data in norm_data_dict.items():
    temp_pka = np.random.normal(6.013, 2.972, 1)[0]
    pka_dict[meta_name] = temp_pka
print(pka_dict)

# save data dict
dest_path_1 = base_path.joinpath("Input/final_1d_data_dict.npy").resolve()
np.save(str(dest_path_1), norm_data_dict)
dest_path_2 = base_path.joinpath("Input/ppm_scale_1d.npy").resolve()
np.save(str(dest_path_2), ppm_scale)
dest_path_3 = base_path.joinpath("Input/pka_dict.npy").resolve()
np.save(str(dest_path_3), pka_dict)


