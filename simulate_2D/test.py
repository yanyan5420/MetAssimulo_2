import pandas as pd
import numpy as np
import argparse
import json
import pathlib
import base64
import io
import re
import plotly.graph_objects as go
from pathlib import Path
import glob
import os
import nmrglue as ng


def read_2d_data(file_path):
    base_path = pathlib.Path(__file__).resolve().parents[1]
    data_path = base_path.joinpath(file_path)
    # print(data_path)

    new_path = Path(data_path)

    dir_list = glob.glob(os.path.join(new_path, '*'))
    # print(dir_list)

    temp_data_dict = dict()
    for sub_dir in dir_list:
        meta_name = str(sub_dir).split("/")[-1]
        dic, data = ng.bruker.read_pdata(sub_dir, shape=(257, 16384))
        temp_data_dict[meta_name] = data
        if meta_name == "Citric Acid":
            print(np.max(data))

    x_scale = np.linspace(13.129, -3.560, 16384)
    y_scale = np.linspace(0.0652, -0.0652, 257)

    return temp_data_dict, x_scale, y_scale


base_path = pathlib.Path(__file__).resolve().parents[1]

cons_df_1 = pd.read_csv(base_path.joinpath("Input/cons_df_1.csv"), index_col=0)
protons_df = pd.read_csv(base_path.joinpath("Input/protons_df.csv"), index_col=0)

with open(base_path.joinpath("Input/hmdb_id_names.json")) as json_file_1:
    hmdb_dict = json.load(json_file_1)



with open(base_path.joinpath("Input/hmdb_normal_concentrations.json")) as json_file_2:
    # pass
    temp = json.load(json_file_2)
    # temp = None


data_dict, x_scale, y_scale = read_2d_data(base_path.joinpath("Input/New_DB_2D/"))
print(len(data_dict))
print(np.max(data_dict["Citric Acid"]))

# print(hex(id(data_dict)))
# print(hex(id(temp)))





print(len(data_dict))
print(len(temp))
print(len(data_dict))

# match_data_dict = db_match_cons(data_dict, cons_df_1, hmdb_dict)
# print(np.max(match_data_dict["citric acid"]))


# with open(base_path.joinpath("Input/hmdb_abnormal_concentrations.json")) as json_file_3:
#     hmdb_abnorm_cons_dict = json.load(json_file_3)
#
# with open(base_path.joinpath("Input/hmdb_id_pka.json")) as json_file_4:
#     hmdb_id_pka_dict = json.load(json_file_4)