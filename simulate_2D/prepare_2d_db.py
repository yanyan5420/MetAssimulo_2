import pandas as pd
import numpy as np
import argparse
import json
import pathlib

from simulate_2D.read_parameters import read_param
from simulate_2D.read_2d_spectra import read_2d_data
from simulate_2D.match_names import db_match_cons, input_match_db
from simulate_2D.preprocess_2d_spectra import remove_water_calibration, filter_noise, smooth_data, normalize_data

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--parameter", dest="parameter_file", required=True, help="Input Parameter Files")
args = parser.parse_args()
param_dict = read_param(args.parameter_file)
file_path_2d = param_dict['file_path_2D']

base_path = pathlib.Path(__file__).resolve().parents[1]
with open(base_path.joinpath("Input/hmdb_id_names.json")) as json_file:
    hmdb_dict = json.load(json_file)
cons_df_1 = pd.read_csv(base_path.joinpath("Input/cons_df_1.csv"), index_col=0)
# cons_df_2 = pd.read_csv("/Users/yanyan/Desktop/Version_0_dash/Input/cons_df_2.csv", index_col=0)
# protons_df = pd.read_csv("/Users/yanyan/Desktop/Version_0_dash/Input/protons_df.csv", index_col=0)

data_dict, x_scale, y_scale = read_2d_data(file_path_2d)  # read 2d data
print(np.max(data_dict["citric acid"]))
match_data_dict = db_match_cons(data_dict, cons_df_1, hmdb_dict)  # format and match db names with cons files

print(np.max(match_data_dict["citric acid"]))

# preprocess 2d data
removed_data_dict = remove_water_calibration(match_data_dict, x_scale, y_scale, [4.5, 5.0], [-0.3, 0.3])
filtered_data_dict = filter_noise(removed_data_dict, 0.05)
smooth_data_dict = smooth_data(filtered_data_dict, 3, 3)
norm_data_dict = normalize_data(smooth_data_dict)
# print(norm_data_dict)
# save data dict
dest_path_1 = base_path.joinpath("Input/final_2d_data_dict.npy").resolve()
np.save(str(dest_path_1), norm_data_dict)
