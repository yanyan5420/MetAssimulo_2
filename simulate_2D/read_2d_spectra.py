import nmrglue as ng
import pandas as pd
import numpy as np
from pathlib import Path
import os
import glob
import pathlib
import re


def read_2d_data(file_path):
    base_path = pathlib.Path(__file__).resolve().parents[1]
    data_path = base_path.joinpath(file_path)

    new_path = Path(data_path)
    dir_list = glob.glob(os.path.join(new_path, '*'))

    data_dict = dict()
    for sub_dir in dir_list:
        meta_name = str(sub_dir).split("/")[-1]
        all_subdirs = [os.path.join(sub_dir, d) for d in os.listdir(sub_dir) if os.path.isdir(os.path.join(sub_dir, d))]
        try:
            dest_dir = [d for d in all_subdirs if d.endswith("1")][0] + "/pdata/1"

        # dest_dir = re.sub(r'/\d+/$', '/1/', sub_dir) + "/pdata/1"
            # print(sub_dir)
            # # print(all_subdirs)
            # print(dest_dir)
            dic, data = ng.bruker.read_pdata(dest_dir, shape=(257, 16384))
            data[-1, :] = data[0, :]
            data_dict[meta_name] = data
        except:
            continue

    x_scale = np.linspace(13.129, -3.560, 16384)
    y_scale = np.linspace(0.0652, -0.0652, 257)

    return data_dict, x_scale, y_scale


def read_2d_cosy(file_path):
    base_path = pathlib.Path(__file__).resolve().parents[1]
    data_path = base_path.joinpath(file_path)

    new_path = Path(data_path)
    dir_list = glob.glob(os.path.join(new_path, '*'))

    data_dict = dict()
    for sub_dir in dir_list:
        meta_name = str(sub_dir).split("/")[-1]
        try:
            dest_dir = sub_dir + "/12/pdata/1"
            dic, data = ng.bruker.read_pdata(dest_dir)
            data_dict[meta_name] = data
        except:
            continue

    x_scale = np.linspace(9.78874, -0.212446, 2048)
    y_scale = np.linspace(9.78874, -0.212446, 2048)

    return data_dict, x_scale, y_scale

